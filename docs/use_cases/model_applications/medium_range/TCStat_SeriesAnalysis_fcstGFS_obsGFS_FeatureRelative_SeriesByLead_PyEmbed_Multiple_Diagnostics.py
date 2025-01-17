"""
Multi_Tool: Feature Relative by Lead using Multiple User-Defined Fields 
========================================================================

model_applications/medium_range/
TCStat_SeriesAnalysis_fcstGFS
_obsGFS_FeatureRelative
_SeriesByLead_PyEmbed_Multiple_Diagnostics.conf

"""

##############################################################################
# Scientific Objective
# --------------------
# This use case calls multiple tools to produce diagnostic plots of systematic erros relative to a
# feature (e.g. hurricane, MCS, etc...). This use case calls two user provided python scripts that
# calculate diagnostics of interest (e.g. integrated vapor transport, potential vorticity, etc...).
# These user diagnostics are then used to define the systematic errors. This example calculates
# statistics over varying forecast leads with the ability to define lead groupings.
# This use case is very similar to the Multi_Tools: Feature Relative by Lead use case and the
# Multi_Tools: Feature Relative by Lead using User-Defined Fields.
# (ADeck,GFS:BDeck,GFS:ATCF,Grib2)
#
# By maintaining focus of each evaluation time (or evaluation time series, in this case)
# on a user-defined area around a cyclone, the model statistical errors associated
# with cyclonic physical features (moisture flux, stability, strength of upper-level
# PV anomaly and jet, etc.) can be related directly to the model forecasts and provide
# improvement guidance by accurately depicting interactions with significant weather
# features around and within the cyclone. This is in contrast to the traditional
# method of regional averaging cyclone observations in a fixed grid, which
# "smooths out" system features and limits the meaningful metrics that can be gathered.
# Specifically, this use case creates bins of forecast lead times as specified by the
# given ranges which provides additional insight directly into forecast lead time accuracy.
#
# Additionally, the ability to calculate model statistical errors based on user provided diagnostics
# allows the user to customize the feature relative analysis to suit their needs.

##############################################################################
# Datasets
# --------
#
# This use case compares the Global Forecast System (GFS) forecast to the GFS analysis for
# hurricane Dorian. It is based on three user provided python scripts that calculate the diagnostic 
# integrated vaport transport (IVT) baroclinic potential vorticity (PV), and saturation equivalent potential temperature (SEPT), respectively. 
# 
#  - Variables required to calculate IVT:
#    Levels required: all pressure levels >= 100mb
#    #. Temperature
#    #. v- component of wind
#    #. u- component of wind
#    #. Geopotential height
#    #. Specific humidity OR Relative Humidity
#
#  - Variables required to calculate PV:
#    Levels required: all pressure levels >= 100mb
#    #. U-wind
#    #. V-wind
#    #. Temperature
#
#  - Variables required to calculate saturation equivalent potential temperature:
#    Levels required: all pressure levels >= 100mb
#    #. Temperature
#
#  - Forecast dataset: GFS Grid 4 Forecast
#    GFS Forecast data can be found at the following website: https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs
#    - Initialization date: 20190830
#    - Initialization hours: 00, 06, 12, 18 UTC
#    - Lead times: 90, 96, 102, 108, 114
#    - Format: Grib2
#    - Resolution: 0.5 degree
#  - Observation dataset: GFS Grid 4 Analysis
#    GFS Analysis data can be found at the following website: https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs
#    - Valid date/time range: 20190902_18 - 20190904_12 every 6 hours
#    - Format: Grib2
#    - Resolution: 0.5 degree
#  - Hurricane Track Data
#    Hurricane track data can be found at the following website: http://hurricanes.ral.ucar.edu/repository/data/
#    - ADeck Track File: aal052019.dat
#    - BDeck Track File: bal052019.dat
#

##############################################################################
# External Dependencies
# ---------------------
#
# You will need to use a version of Python 3.7+ that has the following packages installed:
#
# * netCDF4
# * pygrib
# * cfgrib
# * metpy
# * xarray
#
# If the version of Python used to compile MET did not have these libraries at the time of compilation, you will need to add these packages or create a new Python environment with these packages.
#
# If this is the case, you will need to set the MET_PYTHON_EXE environment variable to the path of the version of Python you want to use. If you want this version of Python to only apply to this use case, set it in the [user_env_vars] section of a METplus configuration file.::
#
#    [user_env_vars]
#    MET_PYTHON_EXE = /path/to/python/with/required/packages/bin/python

##############################################################################
# METplus Components
# ------------------
#
# This use case first runs PyEmbedIngest to run the user provided python scripts to calculate the
# desired diagnostics (in this example, IVT, PV and SEPT). PyEmbedIngest runs the RegridDataPlane tool 
# to write IVT, PV, and SEPTto a MET readable netCDF file. Then TCPairs and ExtractTiles are run to 
# generate matched tropical cyclone data and regrid them into appropriately-sized tiles
# along a storm track. The MET tc-stat tool is used to filter the track data and the MET 
# regrid-dataplane tool is used to regrid the data (GRIB1 or GRIB2 into netCDF). 
# Next, a series analysis by lead time is performed on the results and plots (.ps and .png) are 
# generated for all variable-level-stat combinations from the specified variables, levels, 
# and requested statistics. If lead grouping is turned on, the final results are aggregated into 
# forecast hour groupings as specified by the start, end and increment in the METplus configuration 
# file, as well as labels to identify each forecast hour grouping. If lead grouping is not turned out
# the final results will be written out for each requested lead time.

##############################################################################
# METplus Workflow
# ----------------
#
# This use case loops by process which means that each tool is run for all times before moving to the
# next tool. The tool order is as follows:
# 
# PyEmbedIngest, TCPairs, ExtractTiles, SeriesByLead
#
# This example loops by forecast/lead time (with begin, end, and increment as specified in the METplus
# TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_Multiple_Diagnostics.conf file). 
#
# 4 initialization times will be run over 5 lead times:
#
# | **Init:** 20190830_00Z
# | **Forecast lead:** 90, 96, 102, 108, 114
# |
#
# | **Init:** 20190830_06Z
# | **Forecast lead:** 90, 96, 102, 108, 114
# |
#
# | **Init:** 20190830_12Z
# | **Forecast lead:** 90, 96, 102, 108, 114
# |
#
# | **Init:** 20190830_18Z
# | **Forecast lead:** 90, 96, 102, 108, 114
# |
#

##############################################################################
# METplus Configuration
# ---------------------
#
# METplus first loads all of the configuration files found in parm/metplus_config,
# then it loads any configuration files passed to METplus via the command line
# with the -c option, i.e. -c parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_Multiple_Diagnostics.conf
#
# .. highlight:: bash
# .. literalinclude:: ../../../../parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics.conf
#

#############################################################################
# MET Configuration
# ---------------------
#
# METplus sets environment variables based on user settings in the METplus configuration file. 
# See :ref:`How METplus controls MET config file settings<metplus-control-met>` for more details. 
#
# **YOU SHOULD NOT SET ANY OF THESE ENVIRONMENT VARIABLES YOURSELF! THEY WILL BE OVERWRITTEN BY METPLUS WHEN IT CALLS THE MET TOOLS!**
#
# If there is a setting in the MET configuration file that is currently not supported by METplus you'd like to control, please refer to:
# :ref:`Overriding Unsupported MET config file settings<met-config-overrides>`
#
# **TCPairsConfig_wrapped**
#
# .. note:: See the :ref:`TCPairs MET Configuration<tc-pairs-met-conf>` section of the User's Guide for more information on the environment variables used in the file below:
#
# .. highlight:: bash
# .. literalinclude:: ../../../../parm/met_config/TCPairsConfig_wrapped
#
# **TCStatConfig_wrapped**
#
# .. note:: See the :ref:`TCStat MET Configuration<tc-stat-met-conf>` section of the User's Guide for more information on the environment variables used in the file below:
#
# .. highlight:: bash
# .. literalinclude:: ../../../../parm/met_config/TCStatConfig_wrapped
#
# **SeriesAnalysisConfig_wrapped**
#
# .. note:: See the :ref:`SeriesAnalysis MET Configuration<series-analysis-met-conf>` section of the User's Guide for more information on the environment variables used in the file below:
#
# .. highlight:: bash
# .. literalinclude:: ../../../../parm/met_config/SeriesAnalysisConfig_wrapped

##############################################################################
# Python Embedding
# ----------------
#
# This use case uses four Python embedding scripts to read input data, two for the forecast data and two for the analysis data.
# The multiple datatype input requires the two-script approach.
#
# parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_ivt_fcst.py
#
# .. highlight:: python
# .. literalinclude:: ../../../../parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_ivt_fcst.py
#
# parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_pv_fcst.py
#
# .. highlight:: python
# .. literalinclude:: ../../../../parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_pv_fcst.py
#
# parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_sept_fcst.py
#
# .. highlight:: python
# .. literalinclude:: ../../../../parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_sept_fcst.py
#
#
#
# parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_ivt_analysis.py
#
# .. highlight:: python
# .. literalinclude:: ../../../../parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_ivt_analysis.py
#
# parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_pv_analysis.py
#
# .. highlight:: python
# .. literalinclude:: ../../../../parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_pv_analysis.py
#
# parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_sept_analysis.py
#
# .. highlight:: python
# .. literalinclude:: ../../../../parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics/gfs_sept_analysis.py
#

##############################################################################
# Running METplus
# ---------------
#
# This use case can be run two ways:
#
# 1) Passing in TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics.conf, 
# then a user-specific system configuration file::
#
#        run_metplus.py \
#        /path/to/METplus/parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics.conf \
#        /path/to/user_system.conf
#
# 2) Modifying the configurations in parm/metplus_config, then passing in TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics.conf::
#
#        run_metplus.py \
#        /path/to/METplus/parm/use_cases/model_applications/medium_range/TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multiple_Diagnostics.conf
#
# The former method is recommended. Whether you add them to a user-specific configuration file or modify the metplus_config files, the following variables must be set correctly:
#
# * **INPUT_BASE** - Path to directory where sample data tarballs are unpacked (See Datasets section to obtain tarballs). This is not required to run METplus, but it is required to run the examples in parm/use_cases
# * **OUTPUT_BASE** - Path where METplus output will be written. This must be in a location where you have write permissions
# * **MET_INSTALL_DIR** - Path to location where MET is installed locally
#
#  and for the [exe] section, you will need to define the location of NON-MET executables.
#  If the executable is in the user's path, METplus will find it from the name. 
#  If the executable is not in the path, specify the full path to the executable here (i.e. CONVERT = /usr/bin/convert)
#  The following executables are required for performing series analysis use cases:
#
# Example User Configuration File::
#
#   [dir]
#   INPUT_BASE = /path/to/sample/input/data
#   OUTPUT_BASE = /path/to/output/dir
#   MET_INSTALL_DIR = /path/to/met-X.Y
#
#   [exe]
#   CONVERT = /path/to/convert
#

##############################################################################
# Expected Output
# ---------------
#
# A successful run will output the following both to the screen and to the logfile::
#
#   INFO: METplus has successfully finished running.
#
# Refer to the value set for **OUTPUT_BASE** to find where the output data was generated.
# Output for this use case will be found in subdirectories of the 'series_analysis_lead' directory (relative to **OUTPUT_BASE**):
# 
# * series_animate
# * series_F090
# * series_F096
# * series_F102
# * series_F108
# * series_F114
#
# | The series_animate directory contains the animations of the series analysis in .gif format for all variable, level, and statistics combinations:
#
#    series_animate_<varname>_<level>_<stat>.gif
#
# | The series_FHHH directories contains files that have the following format:
# 
#   ANLY_FILES_FHHH
#
#   FCST_ASCII_FILES_FHHH
#
#   series_FHHH_<varname>_<level>_<stat>.png
#
#   series_FHHH_<varname>_<level>_<stat>.ps
#
#   series_FHHH_<varname>_<level>_<stat>.nc
#
#   Where:
#
#    **HHH** is the forecast hour/lead time in hours
#
#    **varname** is the variable of interest, as specified in the METplus series_by_lead_all_fhrs config file
#
#    **level**  is the level of interest, as specified in the METplus series_by_lead_all_fhrs config file
#
#    **stat** is the statistic of interest, as specified in the METplus series_by_lead_all_fhrs config file.
#

##############################################################################
# Keywords
# --------
#
# .. note::
#
#   * TCPairsToolUseCase
#   * SeriesByLeadUseCase
#   * TCStatToolUseCase
#   * RegridDataPlaneToolUseCase
#   * PyEmbedIngestToolUseCase
#   * MediumRangeAppUseCase
#   * SeriesAnalysisUseCase
#   * GRIB2FileUseCase
#   * FeatureRelativeUseCase
#   * SBUOrgUseCase
#   * DiagnosticsUseCase
#   * RuntimeFreqUseCase
#
#   Navigate to the :ref:`quick-search` page to discover other similar use cases.
#
#
#
# sphinx_gallery_thumbnail_path = '_static/medium_range-TCStat_SeriesAnalysis_fcstGFS_obsGFS_FeatureRelative_SeriesByLead_PyEmbed_Multivariate_Diagnostics.png'
#
