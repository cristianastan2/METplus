"""
Program Name: compare_gridded_wrapper.py
Contact(s): George McCabe
Abstract:
History Log:  Initial version
Usage:
Parameters: None
Input Files:
Output Files:
Condition codes: 0 for success, 1 for failure
"""

import os

from ..util import met_util as util
from ..util import do_string_sub, ti_calculate
from ..util import parse_var_list
from . import CommandBuilder

'''!@namespace CompareGriddedWrapper
@brief Common functionality to wrap similar MET applications
that compare gridded data
Call as follows:
@code{.sh}
Cannot be called directly. Must use child classes.
@endcode
'''

class CompareGriddedWrapper(CommandBuilder):
    """!Common functionality to wrap similar MET applications
that reformat gridded data
    """

    def __init__(self, config, instance=None):
        # set app_name if not set by child class to allow tests to run on this wrapper
        if not hasattr(self, 'app_name'):
            self.app_name = 'compare_gridded'

        super().__init__(config, instance=instance)
        # check to make sure all necessary probabilistic settings are set correctly
        # this relies on the subclass to finish creating the c_dict, so it has to
        # be checked after that happens
        self.check_probabilistic_settings()

    def create_c_dict(self):
        """!Create dictionary from config items to be used in the wrapper
            Allows developer to reference config items without having to know
            the type and consolidates config get calls so it is easier to see
            which config variables are used in the wrapper"""
        c_dict = super().create_c_dict()

        self.add_met_config(name='model',
                            data_type='string',
                            metplus_configs=['MODEL'])

        self.add_met_config(name='obtype',
                            data_type='string',
                            metplus_configs=['OBTYPE'])

        # set old MET config items for backwards compatibility
        c_dict['MODEL_OLD'] = self.config.getraw('config', 'MODEL', 'FCST')
        c_dict['OBTYPE_OLD'] = self.config.getraw('config', 'OBTYPE', 'OBS')

        # INPUT_BASE is not required unless it is referenced in a config file
        # it is used in the use case config files. Don't error if it is not set
        # to a value that contains /path/to
        c_dict['INPUT_BASE'] = self.config.getdir_nocheck('INPUT_BASE', '')

        c_dict['FCST_IS_PROB'] = self.config.getbool('config', 'FCST_IS_PROB', False)
        # if forecast is PROB, get variable to check if prob is in GRIB PDS
        # it can be unset if the INPUT_DATATYPE is NetCDF, so check that after
        # the entire c_dict is created
        if c_dict['FCST_IS_PROB']:
            c_dict['FCST_PROB_IN_GRIB_PDS'] = self.config.getbool('config', 'FCST_PROB_IN_GRIB_PDS', '')

        c_dict['OBS_IS_PROB'] = self.config.getbool('config', 'OBS_IS_PROB', False)
        # see comment for FCST_IS_PROB
        if c_dict['OBS_IS_PROB']:
            c_dict['OBS_PROB_IN_GRIB_PDS'] = self.config.getbool('config', 'OBS_PROB_IN_GRIB_PDS', '')

        c_dict['FCST_PROB_THRESH'] = None
        c_dict['OBS_PROB_THRESH'] = None

        c_dict['ALLOW_MULTIPLE_FILES'] = False
        c_dict['NEIGHBORHOOD_WIDTH'] = ''
        c_dict['NEIGHBORHOOD_SHAPE'] = ''

        self.handle_regrid(c_dict)

        self.handle_description()

        # handle window variables [FCST/OBS]_[FILE_]_WINDOW_[BEGIN/END]
        self.handle_file_window_variables(c_dict)

        self.add_met_config(name='output_prefix',
                            data_type='string')

        c_dict['VAR_LIST_TEMP'] = parse_var_list(self.config,
                                                 met_tool=self.app_name)

        return c_dict

    def set_environment_variables(self, time_info):
        """!Set environment variables that will be read set when running this tool.
            Wrappers could override it to set wrapper-specific values, then call super
            version to handle user configs and printing
            Args:
              @param time_info dictionary containing timing info from current run"""

        self.get_output_prefix(time_info)

        # set old environment variable values for backwards compatibility
        self.add_env_var('MODEL', self.c_dict.get('MODEL_OLD', ''))
        self.add_env_var('OBTYPE', self.c_dict.get('OBTYPE_OLD', ''))
        self.add_env_var('REGRID_TO_GRID',
                         self.c_dict.get('REGRID_TO_GRID', 'NONE'))

        super().set_environment_variables(time_info)

    def check_probabilistic_settings(self):
        """!If dataset is probabilistic, check if *_PROB_IN_GRIB_PDS or INPUT_DATATYPE
            are set. If not enough information is set, report an error and set isOK to False"""
        for dtype in ['FCST', 'OBS']:
            if self.c_dict[f'{dtype}_IS_PROB']:
                # if the data type is NetCDF, then we know how to
                # format the probabilistic fields
                if self.c_dict[f'{dtype}_INPUT_DATATYPE'] != 'GRIB':
                    continue

                # if the data is grib, the user must specify if the data is in
                # the GRIB PDS or not
                if self.c_dict[f'{dtype}_PROB_IN_GRIB_PDS'] == '':
                    self.log_error(f"If {dtype}_IS_PROB is True, you must set {dtype}_PROB_IN"
                                   "_GRIB_PDS unless the forecast datatype is set to NetCDF")
                    self.isOK = False

    def run_at_time(self, input_dict):
        """! Runs the MET application for a given run time. This function loops
              over the list of forecast leads and runs the application for each.
              Args:
                @param input_dict dictionary containing time information
        """

        # loop of forecast leads and process each
        lead_seq = util.get_lead_sequence(self.config, input_dict)
        for lead in lead_seq:
            input_dict['lead'] = lead

            # set current lead time config and environment variables
            time_info = ti_calculate(input_dict)

            self.logger.info("Processing forecast lead {}".format(time_info['lead_string']))

            if util.skip_time(time_info, self.c_dict.get('SKIP_TIMES', {})):
                self.logger.debug('Skipping run time')
                continue

            for custom_string in self.c_dict['CUSTOM_LOOP_LIST']:
                if custom_string:
                    self.logger.info(f"Processing custom string: {custom_string}")

                time_info['custom'] = custom_string

                # Run for given init/valid time and forecast lead combination
                self.run_at_time_once(time_info)

    def run_at_time_once(self, time_info):
        """! Build MET command for a given init/valid time and forecast lead combination
              Args:
                @param time_info dictionary containing timing information
        """
        self.clear()

        var_list = util.sub_var_list(self.c_dict['VAR_LIST_TEMP'],
                                     time_info)

        if not var_list and not self.c_dict.get('VAR_LIST_OPTIONAL', False):
            self.log_error('No input fields were specified. You must set '
                           f'[FCST/OBS]_VAR<n>_[NAME/LEVELS].')
            return None

        if self.c_dict.get('ONCE_PER_FIELD', False):
            # loop over all fields and levels (and probability thresholds) and
            # call the app once for each
            for var_info in var_list:
                self.clear()
                self.c_dict['CURRENT_VAR_INFO'] = var_info
                self.run_at_time_one_field(time_info, var_info)
        else:
            # loop over all variables and all them to the field list, then call the app once
            if var_list:
                self.c_dict['CURRENT_VAR_INFO'] = var_list[0]

            self.run_at_time_all_fields(time_info)

    def run_at_time_one_field(self, time_info, var_info):
        """! Build MET command for a single field for a given
             init/valid time and forecast lead combination
              Args:
                @param time_info dictionary containing timing information
                @param var_info object containing variable information
        """

        # get model to compare, return None if not found
        model_path = self.find_model(time_info,
                                     var_info,
                                     mandatory=True,
                                     return_list=True)
        if model_path is None:
            return

        self.infiles.extend(model_path)
        # get observation to compare, return None if not found
        obs_path, time_info = self.find_obs_offset(time_info,
                                                   var_info,
                                                   mandatory=True,
                                                   return_list=True)
        if obs_path is None:
            return

        self.infiles.extend(obs_path)

        # get field info field a single field to pass to the MET config file
        fcst_field_list = self.get_field_info(v_level=var_info['fcst_level'],
                                              v_thresh=var_info['fcst_thresh'],
                                              v_name=var_info['fcst_name'],
                                              v_extra=var_info['fcst_extra'],
                                              d_type='FCST')

        obs_field_list = self.get_field_info(v_level=var_info['obs_level'],
                                             v_thresh=var_info['obs_thresh'],
                                             v_name=var_info['obs_name'],
                                             v_extra=var_info['obs_extra'],
                                             d_type='OBS')

        if fcst_field_list is None or obs_field_list is None:
            return

        fcst_fields = ','.join(fcst_field_list)
        obs_fields = ','.join(obs_field_list)

        self.format_field('FCST', fcst_fields)
        self.format_field('OBS', obs_fields)

        self.process_fields(time_info)

    def run_at_time_all_fields(self, time_info):
        """! Build MET command for all of the field/level combinations for a given
             init/valid time and forecast lead combination
              Args:
                @param time_info dictionary containing timing information
        """
        var_list = util.sub_var_list(self.c_dict['VAR_LIST_TEMP'],
                                     time_info)

        # get model from first var to compare
        model_path = self.find_model(time_info,
                                     var_list[0],
                                     mandatory=True,
                                     return_list=True)
        if not model_path:
            return

        self.infiles.extend(model_path)
        # get observation to from first var compare
        obs_path, time_info = self.find_obs_offset(time_info,
                                                   var_list[0],
                                                   mandatory=True,
                                                   return_list=True)
        if obs_path is None:
            return

        self.infiles.extend(obs_path)

        fcst_field_list = []
        obs_field_list = []
        for var_info in var_list:
            next_fcst = self.get_field_info(v_level=var_info['fcst_level'],
                                            v_thresh=var_info['fcst_thresh'],
                                            v_name=var_info['fcst_name'],
                                            v_extra=var_info['fcst_extra'],
                                            d_type='FCST')

            next_obs = self.get_field_info(v_level=var_info['obs_level'],
                                           v_thresh=var_info['obs_thresh'],
                                           v_name=var_info['obs_name'],
                                           v_extra=var_info['obs_extra'],
                                           d_type='OBS')

            if next_fcst is None or next_obs is None:
                return

            fcst_field_list.extend(next_fcst)
            obs_field_list.extend(next_obs)

        fcst_field = ','.join(fcst_field_list)
        obs_field = ','.join(obs_field_list)

        self.format_field('FCST', fcst_field)
        self.format_field('OBS', obs_field)

        self.process_fields(time_info)

    def process_fields(self, time_info):
        """! Set and print environment variables, then build/run MET command
              Args:
                @param time_info dictionary with time information
                @param fcst_field field information formatted for MET config file
                @param obs_field field information formatted for MET config file
                @param ens_field field information formatted for MET config file
                only used for ensemble_stat
        """
        # set config file since command is reset after each run
        self.param = do_string_sub(self.c_dict['CONFIG_FILE'],
                                   **time_info)

        self.set_current_field_config()

        # set up output dir with time info
        if not self.find_and_check_output_file(time_info,
                                               is_directory=True):
            return

        # set environment variables needed by MET config file
        self.set_environment_variables(time_info)

        # check if METplus can generate the command successfully
        cmd = self.get_command()
        if cmd is None:
            self.log_error("Could not generate command")
            return

        # run the MET command
        self.build()

    def create_and_set_output_dir(self, time_info):
        """! Builds the full output dir path with valid or init time
              Creates output directory if it doesn't already exist
              Args:
                @param time_info dictionary with time information
        """
        out_dir = self.c_dict['OUTPUT_DIR']

        # use output template if it is set
        # if output template is not set, do not add any extra directories to path
        out_template_name = '{}_OUTPUT_TEMPLATE'.format(self.app_name.upper())
        if self.config.has_option('config',
                                  out_template_name):
            template = self.config.getraw('config',
                                          out_template_name)
            # perform string substitution to get full path
            extra_path = do_string_sub(template,
                                       **time_info)
            out_dir = os.path.join(out_dir, extra_path)

        # create full output dir if it doesn't already exist
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # set output dir for wrapper
        self.outdir = out_dir

    def get_command(self):
        """! Builds the command to run the MET application
           @rtype string
           @return Returns a MET command with arguments that you can run
        """
        if self.app_path is None:
            self.log_error('No app path specified. '
                              'You must use a subclass')
            return None

        cmd = '{} -v {} '.format(self.app_path, self.c_dict['VERBOSITY'])
        for arg in self.args:
            cmd += arg + " "

        if len(self.infiles) == 0:
            self.log_error("No input filenames specified")
            return None

        # add forecast file
        fcst_file = self.infiles[0]
        if fcst_file.startswith('PYTHON'):
            fcst_file = f"'{fcst_file}'"
        cmd += f'{fcst_file} '

        # add observation file
        obs_file = self.infiles[1]
        if obs_file.startswith('PYTHON'):
            obs_file = f"'{obs_file}'"
        cmd += f'{obs_file} '

        if self.param == '':
            self.log_error('Must specify config file to run MET tool')
            return None

        cmd += self.param + ' '

        if self.outdir == "":
            self.log_error("No output directory specified")
            return None

        cmd += '-outdir {}'.format(self.outdir)
        return cmd

    def handle_interp_dict(self, uses_field=False):
        """! Reads config variables for interp dictionary, i.e.
             _INTERP_VLD_THRESH, _INTERP_SHAPE, _INTERP_METHOD, and
             _INTERP_WIDTH. Also _INTERP_FIELD if specified

            @param uses_field if True, read field variable as well
             (default is False)
        """
        items = {
            'vld_thresh': 'float',
            'shape': ('string', 'remove_quotes'),
            'type': ('dict', None, {
                'method': ('string', 'remove_quotes'),
                'width': 'int',
            }),
        }
        if uses_field:
            items['field'] = ('string', 'remove_quotes')

        self.add_met_config_dict('interp', items)
