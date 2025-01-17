"""
RegridDataPlane: Using Python Embedding
=============================================================================

met_tool_wrapper/RegridDataPlane/RegridDataPlane_python
_embedding.conf

"""
##############################################################################
# Scientific Objective
# --------------------
#
# None. Simply regridding data to match a desired grid domain for comparisons.

##############################################################################
# Datasets
# --------
#
# | **Forecast:** ASCII sample file
#
# | **Location:** All of the input data required for this use case can be found in the met_test sample data tarball. Click here to the METplus releases page and download sample data for the appropriate release: https://github.com/dtcenter/METplus/releases
# | This tarball should be unpacked into the directory that you will set the value of INPUT_BASE. See `Running METplus`_ section for more information.
# |

##############################################################################
# METplus Components
# ------------------
#
# This use case utilizes the METplus RegridDataPlane wrapper to generate a command to run the MET tool RegridDataPlane if all required files are found.

##############################################################################
# METplus Workflow
# ----------------
#
# RegridDataPlane is the only tool called in this example. It processes a single run time, but the data does not contain any time information in the filename, so the run time is irrelevant.
#
# This use case regrids data to another domain specified with REGRID_DATA_PLANE_VERIF_GRID. This is done so that
# forecast and observation comparisons are done on the same grid. Many MET comparison tools have regridding capabilities
# built in. However, if the same file is read for comparisons multiple times, it is redundant to regrid that file each time.
# Running RegridDataPlane allows you to regrid once and use the output in many comparisons/evaluations.

##############################################################################
# METplus Configuration
# ---------------------
#
# METplus first loads all of the configuration files found in parm/metplus_config,
# then it loads any configuration files passed to METplus via the command line
# with the -c option, i.e. -c parm/use_cases/met_tool_wrapper/RegridDataPlane/RegridDataPlane_python_embedding.conf
#
# .. highlight:: bash
# .. literalinclude:: ../../../../parm/use_cases/met_tool_wrapper/RegridDataPlane/RegridDataPlane_python_embedding.conf

##############################################################################
# MET Configuration
# -----------------
#
# None. RegridDataPlane does not use configuration files.
#

##############################################################################
# Python Embedding
# ----------------
#
# This use case calls a Python script to read the input data.
# The Python script is stored in the MET repository: /path/to/MET/installation/share/met/python/read_ascii_numpy.py
#
# `read_ascii_numpy.py <https://github.com/dtcenter/MET/blob/develop/met/scripts/python/read_ascii_numpy.py>`_

##############################################################################
# Running METplus
# ---------------
#
# This use case can be run two ways:
#
# 1) Passing in RegridDataPlane_python_embedding.conf then a user-specific system configuration file::
#
#        run_metplus.py -c /path/to/METplus/parm/use_cases/met_tool_wrapper/RegridDataPlane/RegridDataPlane_python_embedding.conf -c /path/to/user_system.conf
#
# 2) Modifying the configurations in parm/metplus_config, then passing in RegridDataPlane_python_embedding.conf::
#
#        run_metplus.py -c /path/to/METplus/parm/use_cases/met_tool_wrapper/RegridDataPlane/RegridDataPlane_python_embedding.conf
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
# Output for this use case will be found in met_tool_wrapper/RegridDataPlane/regrid_py (relative to **OUTPUT_BASE**)
# and will contain the following file:
#
# * numpy_data.nc

##############################################################################
# Keywords
# --------
#
# .. note::
#
#   * RegridDataPlaneToolUseCase
#   * PythonEmbeddingFileUseCase
#
#   Navigate to the :ref:`quick-search` page to discover other similar use cases.
#
#
#
# sphinx_gallery_thumbnail_path = '_static/met_tool_wrapper-RegridDataPlane.png'
#
