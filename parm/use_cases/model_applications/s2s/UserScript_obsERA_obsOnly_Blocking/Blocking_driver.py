#!/usr/bin/env python3
import sys
import os
import numpy as np
import datetime
import netCDF4
import logging
import warnings

# Load the blocking calculation from METcalcpy and plotting from METplotpy
from metcalcpy.contributed.blocking_weather_regime.Blocking import BlockingCalculation
from metcalcpy.contributed.blocking_weather_regime.Blocking_WeatherRegime_util import parse_steps, get_filenames_list, loop_mpr_write, write_mpr_file
from metplotpy.contributed.blocking_s2s import plot_blocking as pb
from metplotpy.contributed.blocking_s2s.CBL_plot import create_cbl_plot


def main():

    # Parse the steps listed in the .conf file to determine which ones to run
    steps_list_fcst,steps_list_obs = parse_steps()

    # If there are no steps listed for the forecast or observation, give a warning that nothing will be run
    if not steps_list_obs and not steps_list_fcst:
        warnings.warn('No processing steps requested for either the model or observations,')
        warnings.warn(' nothing will be run')
        warnings.warn('Set FCST_STEPS and/or OBS_STEPS in the [user_env_vars] section to process data')


    """
    Set up Blocking Calculation, plotting, and mpr output
    """
    # Set up the class for the forecast and observations
    steps_fcst = BlockingCalculation('FCST')
    steps_obs = BlockingCalculation('OBS')

    # Check to see if there is an output plot directory
    # If an output plot directory does not exist, create it as {OUTPUT_BASE}/plots
    oplot_dir = os.environ.get('BLOCKING_PLOT_OUTPUT_DIR','')
    if not oplot_dir:
        obase = os.environ['SCRIPT_OUTPUT_BASE']
        oplot_dir = os.path.join(obase,'plots')
    if not os.path.exists(oplot_dir):
        os.makedirs(oplot_dir)

    # Check to see if there is an output directory for the mpr files
    # If an output mpr directory does not exist, create it as {OUTPUT_BASE}/mpr
    mpr_dir = os.environ.get('BLOCKING_MPR_OUTPUT_DIR','')
    if not mpr_dir:
        obase = os.environ['SCRIPT_OUTPUT_BASE']
        mpr_dir = os.path.join(obase,'mpr')

    # Check to see if CBL's are used from an obs climatology
    use_cbl_obs = os.environ.get('USE_CBL_OBS','False').lower()

    # Get the days per year (season) that should be present
    # The script will fill missing days, but it needs to know how many to expect
    dseasons = int(os.environ['DAYS_PER_SEASON'])

    # Grab the text files that contain listings of the anomaly data to use for the CBL calculation
    obs_cbl_filetxt = os.environ.get('METPLUS_FILELIST_OBS_CBL_INPUT','')
    fcst_cbl_filetxt = os.environ.get('METPLUS_FILELIST_FCST_CBL_INPUT','')

    # Grab the Daily (IBL) text files
    obs_ibl_filetxt = os.environ.get('METPLUS_FILELIST_OBS_IBL_INPUT','')
    fcst_ibl_filetxt = os.environ.get('METPLUS_FILELIST_FCST_IBL_INPUT','')


    """
    Calculate Central Blocking Latitude
    """
    # Check to see if CBL is listed in the observation steps
    if ("CBL" in steps_list_obs):
        logger.logging('Computing Obs CBLs')
        # Get the number of years (seasons) specified for the CBL calculation
        # This is needed to determine the size of the data array for the CBL and IBL calculations
        # The number of years in the CBL calculation can be different from the number of years in the IBL calculation
        cbl_nseasons = int(os.environ['CBL_NUM_SEASONS'])
        # Read in the list of CBL files to use for the observation CBL calculation
        obs_infiles = get_filenames_list(obs_cbl_filetxt)
        # Check to see that the CBL infiles contain the same number of days for all years
        # If not, give an error
        if len(obs_infiles) != (cbl_nseasons*dseasons):
            raise Exception('Invalid Obs data; each year must contain the same date range to calculate seasonal averages.')
        # Run the CBL calculation on the observations
        cbls_obs,lats_obs,lons_obs,mhweight_obs,cbl_time_obs = steps_obs.run_CBL(obs_infiles,cbl_nseasons,dseasons)

    # Check to see if CBL is listed in the forecast steps
    # If the flag to use CBL climatology from the observations is false, calculate CBLs for forecast
    if ("CBL" in steps_list_fcst) and (use_cbl_obs == 'false'):
        # Add in step to use obs for CBLS
        logger.logging('Computing Forecast CBLs')
        # Get the number of years (seasons) specified for the CBL calculation
        # This is needed to determine the size of the data array for the CBL and IBL calculations
        # The number of years in the CBL calculation can be different from the number of years in the IBL calculation
        cbl_nseasons = int(os.environ['CBL_NUM_SEASONS'])
        # Read in the list of CBL files to use for the forecast CBL calculation
        fcst_infiles = get_filenames_list(fcst_cbl_filetxt)
        # Check to see that the CBL infiles contain the same number of days for all years
        # If not, give an error
        if len(fcst_infiles) != (cbl_nseasons*dseasons):
            raise Exception('Invalid Fcst data; each year must contain the same date range to calculate seasonal averages.')
        # Run the CBL calculation on the forecast
        cbls_fcst,lats_fcst,lons_fcst,mhweight_fcst,cbl_time_fcst = steps_fcst.run_CBL(fcst_infiles,cbl_nseasons,dseasons)
        
    # If the flag to use the CBL climatology from the observations is true, then use the obs CBL data for the 
    # forecast CBL data
    elif ("CBL" in steps_list_fcst) and (use_cbl_obs == 'true'):
        # Check to be sure that the observed CBLs were already calculated (listed in the obs_steps)
        if not ("CBL" in steps_list_obs):
            raise Exception('Must run observed CBLs before using them as a forecast.')
        # Use the obs CBL data for the forecast CBL data
        cbls_fcst = cbls_obs
        lats_fcst = lats_obs
        lons_fcst = lons_obs
        mhweight_fcst = mhweight_obs
        cbl_time_fcst = cbl_time_obs


    """
    Plot Central Blocking Latitude
    """
    # Check to see if PLOTCBL is in the observation steps
    if ("PLOTCBL" in steps_list_obs):
        # CBLs need to be calculated before they can be plotted
        # If the CBL step was not in the observation steps, give an error
        if not ("CBL" in steps_list_obs):
            raise Exception('Must run observed CBLs before plotting them.')
        logger.logging('Plotting Obs CBLs')
        # Get the month string and plot output name
        cbl_plot_mthstr = os.environ['OBS_CBL_PLOT_MTHSTR']
        cbl_plot_outname = os.path.join(oplot_dir,os.environ.get('OBS_CBL_PLOT_OUTPUT_NAME','obs_cbl_avg'))
        # Plot the observation CBLs
        create_cbl_plot(lons_obs, lats_obs, cbls_obs, mhweight_obs, cbl_plot_mthstr, cbl_plot_outname, 
            do_averaging=True)
        
    # Check to see if PLOTCBL is the forecast steps
    if ("PLOTCBL" in steps_list_fcst):
        # CBLs need to be calculated before they can be plotted
        # If the CBL step was not in the forecast steps, give an error
        if not ("CBL" in steps_list_fcst):
            raise Exception('Must run forecast CBLs before plotting them.')
        logger.logging('Plotting Forecast CBLs')
        # Get the month string and plot output name
        cbl_plot_mthstr = os.environ['FCST_CBL_PLOT_MTHSTR']
        cbl_plot_outname = os.path.join(oplot_dir,os.environ.get('FCST_CBL_PLOT_OUTPUT_NAME','fcst_cbl_avg'))
        # Plot the forecast CBLs
        create_cbl_plot(lons_fcst, lats_fcst, cbls_fcst, mhweight_fcst, cbl_plot_mthstr, cbl_plot_outname, 
            do_averaging=True)


    """
    Calculate Instantaneously Blocked Longitudes (IBLs)
    """
    # Check to see if IBL is listed in the observation steps
    if ("IBL" in steps_list_obs):
        # Computing IBLs requires the CBL calculation so check to make sure it was in the observation steps
        # If not, give an error
        if not ("CBL" in steps_list_obs):
            raise Exception('Must run observed CBLs before running IBLs.')
        logger.logging('Computing Obs IBLs')
        # Get the number of years (seasons) specified for the IBL calculation
        # This is needed to determine the size of the data array for the CBL and IBL calculations
        # The number of years in the IBL calculation can be different from the number of years in the CBL calculation
        ibl_nseasons = int(os.environ['IBL_NUM_SEASONS'])
        # Read in the list of IBL files to use for the observation IBL calculation
        obs_infiles = get_filenames_list(obs_ibl_filetxt)
        # Check to see that the CBL infiles contain the same number of days for all years
        # If not, give an error
        if len(obs_infiles) != (ibl_nseasons*dseasons):
            raise Exception('Invalid Obs data; each year must contain the same date range to calculate seasonal averages.')
        # Run the IBL calculation on the observation
        ibls_obs,ibl_time_obs = steps_obs.run_Calc_IBL(cbls_obs,obs_infiles,ibl_nseasons,dseasons)
        daynum_obs = np.arange(0,len(ibls_obs[0,:,0]),1)
    
    # Check to see if IBL is listed in the forecast steps
    if ("IBL" in steps_list_fcst):
        # Computing IBLs requries the CBL calculation so check to make sure it was in the forecast steps
        # If it was not, give an error
        if (not "CBL" in steps_list_fcst):
            raise Exception('Must run forecast CBLs or use observed CBLs before running IBLs.')
        logger.logging('Computing Forecast IBLs')
        # Get the number of years (seasons) specified for the IBL calculation
        # This is needed to determine the size of the data array for the CBL and IBL calculations
        # The number of years in the IBL calculation can be different from the number of years in the CBL calculation
        ibl_nseasons = int(os.environ['IBL_NUM_SEASONS'])
        # Read in the list of IBL files to use for the forecast IBL calculation
        fcst_infiles = get_filenames_list(fcst_ibl_filetxt)
        # Check to see that the CBL infiles contain the same number of days for all years
        # If not, give an error
        if len(fcst_infiles) != (ibl_nseasons*dseasons):
            raise Exception('Invalid Fcst data; each year must contain the same date range to calculate seasonal averages.')
        # Run the IBL calculation on the forecast
        ibls_fcst,ibl_time_fcst = steps_fcst.run_Calc_IBL(cbls_fcst,fcst_infiles,ibl_nseasons,dseasons)
        daynum_fcst = np.arange(0,len(ibls_fcst[0,:,0]),1)

    """
    Write out matched pair files for the IBLs
    if IBLs are calculated in both the obs and forecast
    """
    if ("IBL" in steps_list_obs) and ("IBL" in steps_list_fcst):
        # Make sure an output directory is created
        i_mpr_outdir = os.path.join(mpr_dir,'IBL')
        if not os.path.exists(i_mpr_outdir):
            os.makedirs(i_mpr_outdir)
        modname = os.environ.get('MODEL_NAME','GFS')
        maskname = os.environ.get('MASK_NAME','FULL')
        # Create output file prefix
        ibl_outfile_prefix = os.path.join(i_mpr_outdir,'IBL_stat_'+modname)
        # Average cbls for use as latitude
        cbls_avg = np.nanmean(cbls_obs,axis=0)
        # Write Output File
        loop_mpr_write(ibls_obs,ibls_fcst,cbls_avg,lons_obs,ibl_time_obs,ibl_time_fcst,modname,
            'NA','IBLs','block','Z500','IBLs','block','Z500',maskname,'500',ibl_outfile_prefix)


    """
    Plot IBLs
    """
    if("PLOTIBL" in steps_list_obs) and not ("PLOTIBL" in steps_list_fcst):
        # Check to be sure the observation IBL step has been run
        if not ("IBL" in steps_list_obs):
            raise Exception('Must run observed IBLs before plotting them.')
        logger.logging('Plotting Obs IBLs')
        # Get title, output name, and plot label
        ibl_plot_title = os.environ.get('OBS_IBL_PLOT_TITLE','Instantaneous Blocked Longitude')
        ibl_plot_outname = os.path.join(oplot_dir,os.environ.get('OBS_IBL_PLOT_OUTPUT_NAME','obs_IBL_Freq'))
        ibl_plot_label1 = os.environ.get('IBL_PLOT_OBS_LABEL','')
        # Plot observed IBLs
        pb.plot_ibls(ibls_obs,lons_obs,ibl_plot_title,ibl_plot_outname,label1=ibl_plot_label1)
    elif ("PLOTIBL" in steps_list_fcst) and not ("PLOTIBL" in steps_list_obs):
        # Check to be sure the forecast IBL step has ben run
        if not ("IBL" in steps_list_fcst):
            raise Exception('Must run forecast IBLs before plotting them.')
        logger.logging('Plotting Forecast IBLs')
        # Get title, output name, and plot label
        ibl_plot_title = os.environ.get('FCST_IBL_PLOT_TITLE','Instantaneous Blocked Longitude')
        ibl_plot_outname = os.path.join(oplot_dir,os.environ.get('FCST_IBL_PLOT_OUTPUT_NAME','fcst_IBL_Freq'))
        ibl_plot_label1 = os.environ.get('IBL_PLOT_FCST_LABEL','')
        # Plot forecast IBLs
        pb.plot_ibls(ibls_fcst,lons_fcst,ibl_plot_title,ibl_plot_outname,label1=ibl_plot_label1)
    # If PLOTIBL is a step in both the observation and forecast, then put them on the same plot
    elif ("PLOTIBL" in steps_list_obs) and ("PLOTIBL" in steps_list_fcst):
        # Check to be sure the observation and forecast IBL steps have been run
        if (not "IBL" in steps_list_obs) and (not "IBL" in steps_list_fcst):
            raise Exception('Must run forecast and observed IBLs before plotting them.')
        logger.logging('Plotting Obs and Forecast IBLs')
        # Get title, output name, and plot labels
        ibl_plot_title = os.environ['IBL_PLOT_TITLE']
        ibl_plot_outname = os.path.join(oplot_dir,os.environ.get('IBL_PLOT_OUTPUT_NAME','IBL_Freq'))
        # Check to see if there are plot legend labels
        ibl_plot_label1 = os.environ.get('IBL_PLOT_OBS_LABEL','Observation')
        ibl_plot_label2 = os.environ.get('IBL_PLOT_FCST_LABEL','Forecast')
        # Plot observed and forecast IBLs
        pb.plot_ibls(ibls_obs,lons_obs,ibl_plot_title,ibl_plot_outname,data2=ibls_fcst,lon2=lons_fcst,
            label1=ibl_plot_label1,label2=ibl_plot_label2)


    """
    Caluclate Group Instantaneously Blocked Longitudes (GIBLs)
    """
    if ("GIBL" in steps_list_obs):
        # Check to be sure the observed IBL step was run
        if not ("IBL" in steps_list_obs):
            raise Exception('Must run observed IBLs before running GIBLs.')
        logger.logging('Computing Obs GIBLs')
        # Calculate Observed GIBLs
        gibls_obs = steps_obs.run_Calc_GIBL(ibls_obs,lons_obs)

    if ("GIBL" in steps_list_fcst):
        # Check to be sure the forecast IBL step was run
        if not ("IBL" in steps_list_fcst):
            raise Exception('Must run Forecast IBLs before running GIBLs.')
        logger.logging('Computing Forecast GIBLs')
        # Calculate Forecast GIBLs
        gibls_fcst = steps_fcst.run_Calc_GIBL(ibls_fcst,lons_fcst)


    """
    Calculate Blocks
    """
    if ("CALCBLOCKS" in steps_list_obs):
        # Check to be sure the observed GIBL step was run
        if not ("GIBL" in steps_list_obs):
            raise Exception('Must run observed GIBLs before calculating blocks.')
        logger.logging('Computing Obs Blocks')
        # Calculate observed blocks
        block_freq_obs = steps_obs.run_Calc_Blocks(ibls_obs,gibls_obs,lons_obs,daynum_obs)

    if ("CALCBLOCKS" in steps_list_fcst):
        # Check to be sure the forecast GIBL step was run
        if not ("GIBL" in steps_list_fcst):
            raise Exception('Must run Forecast GIBLs before calculating blocks.')
        logger.logging('Computing Forecast Blocks')
        # Calculate forecast blocks
        block_freq_fcst = steps_fcst.run_Calc_Blocks(ibls_fcst,gibls_fcst,lons_fcst,daynum_fcst)


    """ 
    Write out matched pair files for the Blocking
    if blocks are calculated in both the observation and forecast
    """
    if ("CALCBLOCKS" in steps_list_obs) and ("CALCBLOCKS" in steps_list_fcst):
        # Create and Output directory if it doesn't exist
        b_mpr_outdir = os.path.join(mpr_dir,'Blocks')
        if not os.path.exists(b_mpr_outdir):
            os.makedirs(b_mpr_outdir)
        # Get model name and mask name
        modname = os.environ.get('MODEL_NAME','GFS')
        maskname = os.environ.get('MASK_NAME','FULL')
        # Create output file prefix
        blocks_outfile_prefix = os.path.join(b_mpr_outdir,'blocking_stat_'+modname)
        # Average cbls for use as latitude
        cbls_avg = np.nanmean(cbls_obs,axis=0)
        # Write output mpr file
        loop_mpr_write(block_freq_obs,block_freq_fcst,cbls_avg,lons_obs,ibl_time_obs,ibl_time_fcst,modname,
            'NA','Blocks','block','Z500','Blocks','block','Z500',maskname,'500',blocks_outfile_prefix)


    """
    Plot Blocking Frequency
    """
    if ("PLOTBLOCKS" in steps_list_obs):
        # Check to be sure the observed blocking step was run
        if not ("CALCBLOCKS" in steps_list_obs):
            raise Exception('Must compute observed blocks before plotting them.')
        logger.logging('Plotting Obs Blocks')
        # Get plot title and output name
        blocking_plot_title = os.environ.get('OBS_BLOCKING_PLOT_TITLE','Obs Blocking Frequency')
        blocking_plot_outname = os.path.join(oplot_dir,os.environ.get('OBS_BLOCKING_PLOT_OUTPUT_NAME','obs_Block_Freq'))
        # plot observed blocking frequency
        pb.plot_blocks(block_freq_obs,gibls_obs,ibls_obs,lons_obs,blocking_plot_title,blocking_plot_outname)
    if ("PLOTBLOCKS" in steps_list_fcst):
        # Check to be sure the forecast blocking step was run
        if not ("CALCBLOCKS" in steps_list_fcst):
            raise Exception('Must compute forecast blocks before plotting them.')
        logger.logging('Plotting Forecast Blocks')
        # Get plot title and output name
        blocking_plot_title = os.environ.get('FCST_BLOCKING_PLOT_TITLE','Forecast Blocking Frequency')
        blocking_plot_outname = os.path.join(oplot_dir,os.environ.get('FCST_BLOCKING_PLOT_OUTPUT_NAME','fcst_Block_Freq'))
        # plot forecast blocking frequency
        pb.plot_blocks(block_freq_fcst,gibls_fcst,ibls_fcst,lons_fcst,blocking_plot_title,blocking_plot_outname)


if __name__ == "__main__":
    main()
