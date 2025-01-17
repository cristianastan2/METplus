"""
Grid-Stat: Surrogate Severe and Practically Perfect Probabilistic Evaluation 
============================================================================

model_applications/
convection_allowing_models/
GridStat_fcstHRRR_obsPracPerfect
_SurrogateSevereProb.conf

"""

##############################################################################
# Scientific Objective
# --------------------
#
# To evaluate the surrogate severe forecasts at predicting Severe weather
# using the (12Z - 12Z) practically perfect storm reports an obtain 
# probabilistic output statistics.

##############################################################################
# Datasets
# --------
#
#  * Forecast dataset: HRRR Surrogate Severe Data
#  * Observation dataset: Practically Perfect from Local Storm Reports
#

##############################################################################
# METplus Components
# ------------------
#
# This use case runs grid_stat to create probabilistic statistics on 
# surrogate severe from the HRRR model and Practially Perfect observations 
# computed from local storm reports.  

##############################################################################
# METplus Workflow
# ----------------
#
# The grid_stat tool is run for each time. This example loops by valid time.  It
# processes 1 valid time, listed below.
#
# | **Valid:** 2020-02-06_12Z
# | **Forecast lead:** 36
# |

##############################################################################
# METplus Configuration
# ---------------------
#
# METplus first loads all of the configuration files found in parm/metplus_config,
# then it loads any configuration files passed to METplus via the command line
# with the -c option, i.e. -c parm/use_cases/model_applications/convection_allowing_models/GridStat_fcstHRRR_obsPracPerfect_SurrogateSevere.conf
#
# .. highlight:: bash
# .. literalinclude:: ../../../../parm/use_cases/model_applications/convection_allowing_models/GridStat_fcstHRRR_obsPracPerfect_SurrogateSevereProb.conf

##############################################################################
# MET Configuration
# -----------------
#
# METplus sets environment variables based on user settings in the METplus configuration file. 
# See :ref:`How METplus controls MET config file settings<metplus-control-met>` for more details. 
#
# **YOU SHOULD NOT SET ANY OF THESE ENVIRONMENT VARIABLES YOURSELF! THEY WILL BE OVERWRITTEN BY METPLUS WHEN IT CALLS THE MET TOOLS!**
#
# If there is a setting in the MET configuration file that is currently not supported by METplus you'd like to control, please refer to:
# :ref:`Overriding Unsupported MET config file settings<met-config-overrides>`
#
# .. note:: See the :ref:`GridStat MET Configuration<grid-stat-met-conf>` section of the User's Guide for more information on the environment variables used in the file below:
#
# .. highlight:: bash
# .. literalinclude:: ../../../../parm/met_config/GridStatConfig_wrapped

##############################################################################
# Running METplus
# ---------------
#
# This use case can be run two ways:
#
# 1) Passing in GridStat_fcstHRRR_obsPracPerfect_SurrogateSevere.conf then a user-specific system configuration file::
#
#        run_metplus.py -c /path/to/METplus/parm/use_cases/model_applications/convection_allowing_models/GridStat_fcstHRRR_obsPracPerfect_SurrogateSevereProb.conf -c /path/to/user_system.conf
#
# 2) Modifying the configurations in parm/metplus_config, then passing in GridStat_fcstHRRR_obsPracPerfect_SurrogateSevere.conf::
#
#        run_metplus.py -c /path/to/METplus/parm/use_cases/model_applications/convection_allowing_models/GridStat_fcstHRRR_obsPracPerfect_SurrogateSevereProb.conf
#
# The former method is recommended. Whether you add them to a user-specific configuration file or modify the metplus_config files, the following variables must be set correctly:
#
# * **INPUT_BASE** - Path to directory where sample data tarballs are unpacked (See Datasets section to obtain tarballs). This is not required to run METplus, but it is required to run the examples in parm/use_cases
# * **OUTPUT_BASE** - Path where METplus output will be written. This must be in a location where you have write permissions
# * **MET_INSTALL_DIR** - Path to location where MET is installed locally
#
# Example User Configuration File::
#
#   [dir]
#   INPUT_BASE = /path/to/sample/input/data
#   OUTPUT_BASE = /path/to/output/dir
#   MET_INSTALL_DIR = /path/to/met-X.Y 
#
# **NOTE:** All of these items must be found under the [dir] section.

##############################################################################
# Expected Output
# ---------------
#
# A successful run will output the following both to the screen and to the logfile::
#
#   INFO: METplus has successfully finished running.
#
# Refer to the value set for **OUTPUT_BASE** to find where the output data was generated.
# Output for this use case will be found in model_applications/convection_allowing_models/surrogate_severe_prac_perfect/grid_stat/prob (relative to **OUTPUT_BASE**)
# and will contain the following files:
#
# grid_stat_360000L_20200206_120000V_pct.txt
# grid_stat_360000L_20200206_120000V_pjc.txt
# grid_stat_360000L_20200206_120000V_prc.txt
# grid_stat_360000L_20200206_120000V_pstd.txt
# grid_stat_360000L_20200206_120000V.stat

##############################################################################
# Keywords
# --------
#
#
# .. note::
#
#   * GridStatToolUseCase
#   * ConvectionAllowingModelsAppUseCase  
#   * NetCDFFileUseCase 
#   * NOAAHWTOrgUseCase 
#   * NCAROrgUseCase 
#   * NOAAHMTOrgUseCase 
#
#   Navigate to the :ref:`quick-search` page to discover other similar use cases.
#
#
#
# sphinx_gallery_thumbnail_path = '_static/convection_allowing_models-SS_PP_prob.png'
#
