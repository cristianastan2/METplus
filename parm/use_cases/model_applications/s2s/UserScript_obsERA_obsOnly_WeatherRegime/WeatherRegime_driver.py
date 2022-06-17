#!/usr/bin/env python3
import sys
import os
import numpy as np
import netCDF4
import warnings
import logging

# Load the blocking calculation from METcalcpy and plotting from METplotpy
from metcalcpy.contributed.blocking_weather_regime.WeatherRegime import WeatherRegimeCalculation, reorder_fcst_regimes,reorder_fcst_regimes_correlate
from metcalcpy.contributed.blocking_weather_regime.Blocking_WeatherRegime_util import parse_steps, get_filenames_list, read_nc_met, loop_mpr_write, write_mpr_file
from metplotpy.contributed.weather_regime import plot_weather_regime as pwr


"""
Function that reads in the list of data files, and then reads in the data
"""
def read_wr_data(intype,data_filetext,nseasons,dseasons):
    data_infiles = get_filenames_list(data_filetext)
    # Check to be sure the same number of days is given for each year
    if len(data_infiles) != (nseasons*dseasons):
            raise Exception('Invalid '+intype+' data; each year must contain the same date range to calculate seasonal averages.')
    data_invar = os.environ.get(intype+'_WR_VAR','')
    # Read in the data
    z500_in,lats_in,lons_in,timedict_in = read_nc_met(data_infiles,data_invar,nseasons,dseasons)

    return z500_in,lats_in,lons_in,timedict_in


"""
Blocking Driver that calls requested steps
The code for the requested steps is in METcalcpy and METplotpy
"""
def main():

    # Parse the steps listed in the .conf file to determine which ones to run
    steps_list_fcst,steps_list_obs = parse_steps()

    # If there are no steps listed for the forecast or observation, give a warning that nothing will be run
    if not steps_list_obs and not steps_list_fcst:
        warnings.warn('No processing steps requested for either the model or observations,')
        warnings.warn(' nothing will be run')
        warnings.warn('Set FCST_STEPS and/or OBS_STEPS in the [user_env_vars] section to process data')


    """
    Set up Weather Regime Calculation, Plotting, and mpr output
    """
    # Set up the class for the forecast and observations
    steps_obs = WeatherRegimeCalculation('OBS')
    steps_fcst = WeatherRegimeCalculation('FCST')

    # Check to see if there is an output plot directory
    # If an output plot directory does not exist, create it as {OUTPUT_BASE}/plots
    oplot_dir = os.environ.get('WR_PLOT_OUTPUT_DIR','')
    obase = os.environ['SCRIPT_OUTPUT_BASE']
    if not oplot_dir:
        oplot_dir = os.path.join(obase,'plots')
    if not os.path.exists(oplot_dir):
        os.makedirs(oplot_dir)

    # Check to see if there is an output directory for the mpr files
    # If an output mpr directory does not exist, create it as {OUTPUT_BASE}/mpr
    mpr_outdir = os.environ.get('WR_MPR_OUTPUT_DIR','')
    if not mpr_outdir:
        mpr_outdir = os.path.join(obase,'mpr')

    # Get the days per year (season) and number of years (seasons) that should be present
    # The script will fill missing days, but it needs to know how many to expect
    nseasons = int(os.environ['NUM_SEASONS'])
    dseasons = int(os.environ['DAYS_PER_SEASON'])


    """
    Read in the observation and forecast data
    """
    # If observation steps are requested, read in the observation data and detrend
    if ("ELBOW" in steps_list_obs) or ("EOF" in steps_list_obs) or ("KMEANS" in steps_list_obs):
        logging.info('Reading in observation data')
        # Grab the text file that contains listings of the daily averaged data to use for the calculation
        obs_wr_filetxt = os.environ.get('METPLUS_FILELIST_OBS_INPUT','')
        # Read in the data and detrend
        z500_obs,lats_obs,lons_obs,timedict_obs = read_wr_data('OBS',obs_wr_filetxt,nseasons,dseasons)
        z500_detrend_obs,z500_detrend_2d_obs = steps_obs.weights_detrend(lats_obs,lons_obs,z500_obs)

    # If forecast steps are requested, read in the observation data and detrend
    if ("ELBOW" in steps_list_fcst) or ("EOF" in steps_list_fcst) or("KMEANS" in steps_list_fcst):
        logging.info('Reading in Forecast data')
        # Grab the text file that contains listings of the daily averaged data to use for the calculation
        fcst_wr_filetxt = os.environ.get('METPLUS_FILELIST_FCST_INPUT','')
        # Read in the data and detrend
        z500_fcst,lats_fcst,lons_fcst,timedict_fcst = read_wr_data('FCST',fcst_wr_filetxt,nseasons,dseasons)
        z500_detrend_fcst,z500_detrend_2d_fcst = steps_fcst.weights_detrend(lats_fcst,lons_fcst,z500_fcst)


    """
    Compute bend in the elbow to get the optimal number of clusters
    """
    if ("ELBOW" in steps_list_obs):
        logging.info('Running Obs Elbow')
        # Compute the observation elbow
        K_obs,d_obs,mi_obs,line_obs,curve_obs = steps_obs.run_elbow(z500_detrend_2d_obs)

    if ("ELBOW" in steps_list_fcst):
        logging.info('Running Forecast Elbow')
        # Compute the forecast elbow
        K_fcst,d_fcst,mi_fcst,line_fcst,curve_fcst = steps_fcst.run_elbow(z500_detrend_2d_fcst)


    """
    Plot the Elbow
    """
    if ("PLOTELBOW" in steps_list_obs):
        # Check to be sure the observation ELBOW step has been run
        if not ("ELBOW" in steps_list_obs):
            raise Exception('Must run observed Elbow before plotting observed elbow.')
        logging.info('Creating Obs Elbow plot')
        # Get plot title, output name, and plot
        elbow_plot_title = os.environ.get('OBS_ELBOW_PLOT_TITLE','Elbow Method For Optimal k')
        elbow_plot_outname = os.path.join(oplot_dir,os.environ.get('OBS_ELBOW_PLOT_OUTPUT_NAME','obs_elbow'))
        pwr.plot_elbow(K_obs,d_obs,mi_obs,line_obs,curve_obs,elbow_plot_title,elbow_plot_outname)

    if ("PLOTELBOW" in steps_list_fcst):
        # Check to be sure the forecast ELBOW step has been run
        if not ("ELBOW" in steps_list_fcst):
            raise Exception('Must run forecast Elbow before plotting forecast elbow.')
        logging.info('Creating Forecast Elbow plot')
        # Get plot title, output name, and plot
        elbow_plot_title = os.environ.get('FCST_ELBOW_PLOT_TITLE','Elbow Method For Optimal k')
        elbow_plot_outname = os.path.join(oplot_dir,os.environ.get('FCST_ELBOW_PLOT_OUTPUT_NAME','fcst_elbow'))
        pwr.plot_elbow(K_fcst,d_fcst,mi_fcst,line_fcst,curve_fcst,elbow_plot_title,elbow_plot_outname)


    """
    Compute EOFs and Reconstruct heights
    """
    if ("EOF" in steps_list_obs):
        logging.info('Running Obs EOF')
        # Compute observation EOFs and reconstruct heights
        eof_obs,pc_obs,wrnum_obs,variance_fractions_obs = steps_obs.Calc_EOF(z500_obs)
        z500_detrend_2d_obs = steps_obs.reconstruct_heights(eof_obs,pc_obs,z500_detrend_2d_obs.shape)

    if ("EOF" in steps_list_fcst):
        logging.info('Running Forecast EOF')
        # Compute forecast EOFs and reconstruct heights
        eof_fcst,pc_fcst,wrnum_fcst,variance_fractions_fcst = steps_fcst.Calc_EOF(z500_fcst)
        z500_detrend_2d_fcst = steps_fcst.reconstruct_heights(eof_fcst,pc_fcst,z500_detrend_2d_fcst.shape)


    """
    Plot EOFs
    """
    if ("PLOTEOF" in steps_list_obs):
        # Check to be sure the observation EOF step has been run
        if not ("EOF" in steps_list_obs):
            raise Exception('Must run observed EOFs before plotting observed EOFs.')
        logging.info('Plotting Obs EOFs')
        # Get plot levels, output name, and plot
        pltlvls_str = os.environ['EOF_PLOT_LEVELS'].split(',')
        pltlvls = [float(pp) for pp in pltlvls_str]
        eof_plot_outname = os.path.join(oplot_dir,os.environ.get('OBS_EOF_PLOT_OUTPUT_NAME','obs_eof'))
        pwr.plot_eof(eof_obs,wrnum_obs,variance_fractions_obs,lons_obs,lats_obs,eof_plot_outname,pltlvls)

    if ("PLOTEOF" in steps_list_fcst):
        # Check to be sure the forecast EOF step has been run
        if not ("EOF" in steps_list_fcst):
            raise Exception('Must run forecast EOFs before plotting forecast EOFs.')
        logging.info('Plotting Forecast EOFs')
        # Get plot levels, output name, and plot
        pltlvls_str = os.environ['EOF_PLOT_LEVELS'].split(',')
        pltlvls = [float(pp) for pp in pltlvls_str]
        eof_plot_outname = os.path.join(oplot_dir,os.environ.get('FCST_EOF_PLOT_OUTPUT_NAME','fcst_eof'))
        pwr.plot_eof(eof_fcst,wrnum_fcst,variance_fractions_fcst,lons_fcst,lats_fcst,eof_plot_outname,pltlvls)


    """
    Compute K-means and write a matched pair file if both obs and fcst
    K-means are run
    """
    # Check to see if KMEANS is in the observation steps
    if ("KMEANS" in steps_list_obs):
        logging.info('Running Obs K Means')
        # Calcuate K-means and write output file with weather regime classification for each day
        kmeans_obs,wrnum_obs,perc_obs,wrc_obs= steps_obs.run_K_means(z500_detrend_2d_obs,timedict_obs,z500_obs.shape)
        steps_obs.write_K_means_file(timedict_obs,wrc_obs)

    # Check to see if KMEANS is in the forecast steps
    if ("KMEANS" in steps_list_fcst):
        logging.info('Running Forecast K Means')
        # Calcuate K-means
        kmeans_fcst,wrnum_fcst,perc_fcst,wrc_fcst = steps_fcst.run_K_means(z500_detrend_2d_fcst,timedict_fcst,
            z500_fcst.shape)
        # Reorder the forecast regime numbers so the numbers match the observation
        reorder_fcst = os.environ.get('REORDER_FCST','False').lower()
        reorder_fcst_manual = os.environ.get('REORDER_FCST_MANUAL','False').lower()
        # Use spatial correlation to reorder if KMEANS was computed in the observation step
        if (reorder_fcst == 'true') and ("KMEANS" in steps_list_obs): 
            kmeans_fcst,perc_fcst,wrc_fcst = reorder_fcst_regimes_correlate(kmeans_obs,kmeans_fcst,perc_fcst,wrc_fcst,wrnum_fcst)
        # Reorder manually if requested
        if reorder_fcst_manual == 'true':
            fcst_order_str = os.environ['FCST_ORDER'].split(',')
            fcst_order = [int(fo) for fo in fcst_order_str]
            kmeans_fcst,perc_fcst,wrc_fcst = reorder_fcst_regimes(kmeans_fcst,perc_fcst,wrc_fcst,wrnum_fcst,fcst_order)
        # Write output file with weather regime classification for each day
        steps_fcst.write_K_means_file(timedict_fcst,wrc_fcst)

    # Write and output matched pair file if KMEANS was computed for the forecast and observation
    if ("KMEANS" in steps_list_obs) and ("KMEANS" in steps_list_fcst):
        # Get input model, mask name
        modname = os.environ.get('MODEL_NAME','GFS')
        maskname = os.environ.get('MASK_NAME','FULL')
        wrc_obs_mpr = wrc_obs[:,:,np.newaxis]
        wrc_fcst_mpr = wrc_fcst[:,:,np.newaxis]
        # Create output directory if it doesn't exist
        mpr_full_outdir = os.path.join(mpr_outdir,'WeatherRegime')
        if not os.path.exists(mpr_full_outdir):
            os.makedirs(mpr_full_outdir)
        # Write mpr file
        wr_outfile_prefix = os.path.join(mpr_full_outdir,'weather_regime_stat_'+modname)
        loop_mpr_write(wrc_obs_mpr,wrc_fcst_mpr,[0.0],[0.0],timedict_obs,timedict_fcst,modname,'NA',
            'WeatherRegimeClass','class','Z500','WeatherRegimeClass','class','Z500',maskname,'500',
            wr_outfile_prefix)


    """
    Plot KMEANS Weather Regime patterns
    """
    if ("PLOTKMEANS" in steps_list_obs):
        # Check to be sure the observation KMEANS step has been run
        if not ("KMEANS" in steps_list_obs):
            raise Exception('Must run observed Kmeans before plotting observed Kmeans.')
        logging.info('Plotting Obs K Means')
        # Get plot levels, output name, and plot
        pltlvls_str = os.environ['KMEANS_PLOT_LEVELS'].split(',')
        pltlvls = [float(pp) for pp in pltlvls_str]
        kmeans_plot_outname = os.path.join(oplot_dir,os.environ.get('OBS_KMEANS_PLOT_OUTPUT_NAME','obs_kmeans'))
        pwr.plot_K_means(kmeans_obs,wrnum_obs,lons_obs,lats_obs,perc_obs,kmeans_plot_outname,pltlvls)

    if ("PLOTKMEANS" in steps_list_fcst):
        # Check to be sure the forecast KMEANS step has been run
        if not ("KMEANS" in steps_list_fcst):
            raise Exception('Must run forecast Kmeans before plotting forecast Kmeans.')
        logging.info('Plotting Forecast K Means')
        # Get plot levels, output name, and plot
        pltlvls_str = os.environ['KMEANS_PLOT_LEVELS'].split(',')
        pltlvls = [float(pp) for pp in pltlvls_str]
        kmeans_plot_outname = os.path.join(oplot_dir,os.environ.get('FCST_KMEANS_PLOT_OUTPUT_NAME','fcst_kmeans'))
        pwr.plot_K_means(kmeans_fcst,wrnum_fcst,lons_fcst,lats_fcst,perc_fcst,kmeans_plot_outname,pltlvls)


    """
    Calculate Time Frequency and write a matched pair file if both obs
    and fcst Time Frequency are run
    """
    if ("TIMEFREQ" in steps_list_obs):
        # Check to be sure the observation KMEANS step has been run
        if not ("KMEANS" in steps_list_obs):
            raise Exception('Must run observed Kmeans before running frequencies.')
        # Compute the time frequency
        wrfreq_obs,dlen_obs,ts_diff_obs = steps_obs.compute_wr_freq(wrc_obs)

    if ("TIMEFREQ" in steps_list_fcst):
        # Check to be sure the forecast KMEANS step has been run
        if not ("KMEANS" in steps_list_fcst):
            raise Exception('Must run forecast Kmeans before running frequencies.')
        # Compute the time frequency
        wrfreq_fcst,dlen_fcst,ts_diff_fcst = steps_fcst.compute_wr_freq(wrc_fcst)

    # Write and output matched pair file if TIMEFREQ was computed for the forecast and observation
    if ("TIMEFREQ" in steps_list_obs) and ("TIMEFREQ" in steps_list_fcst):
        # Get model name and mask name
        modname = os.environ.get('MODEL_NAME','GFS')
        maskname = os.environ.get('MASK_NAME','FULL')
        wrfreq_obs_mpr = wrfreq_obs[:,:,:,np.newaxis]
        wrfreq_fcst_mpr = wrfreq_fcst[:,:,:,np.newaxis]
        # Set up a time dictionary
        timedict_obs_mpr = {'init':timedict_obs['init'][:,ts_diff_obs-1:],
            'valid':timedict_obs['valid'][:,ts_diff_obs-1:],'lead':timedict_obs['lead'][:,ts_diff_obs-1:]}
        timedict_fcst_mpr = {'init':timedict_fcst['init'][:,ts_diff_fcst-1:],
            'valid':timedict_fcst['valid'][:,ts_diff_fcst-1:],'lead':timedict_fcst['lead'][:,ts_diff_fcst-1:]}
        # Create output directory if it doesn't exist
        mpr_full_outdir = os.path.join(mpr_outdir,'freq')
        if not os.path.exists(mpr_full_outdir):
            os.makedirs(mpr_full_outdir)
        # Write output mpr file
        for wrn in np.arange(wrnum_obs):
            wr_outfile_prefix = os.path.join(mpr_full_outdir,'weather_regime'+str(wrn+1).zfill(2)+'_freq_stat_'+modname)
            loop_mpr_write(wrfreq_obs_mpr[wrn,:,:,:],wrfreq_fcst_mpr[wrn,:,:,:],[0.0],[0.0],timedict_obs,
                timedict_fcst,modname,str(wrn+1).zfill(2),'WeatherRegimeFreq','percent','Z500','WeatherRegimeFreq',
                'percent','Z500',maskname,'500',wr_outfile_prefix)


    """
    Plot Time Frequency
    """
    if ("PLOTFREQ" in steps_list_obs):
        # Check to be sure the observation TIMEFREQ step has been run
        if not ("TIMEFREQ" in steps_list_obs):
            raise Exception('Must run observed Frequency calculation before plotting the frequencies.')
        # Get plot title, output name, compute mean, and plot
        freq_plot_title = os.environ.get('OBS_FREQ_PLOT_TITLE','Seasonal Cycle of WR Days/Week')
        freq_plot_outname = os.path.join(oplot_dir,os.environ.get('OBS_FREQ_PLOT_OUTPUT_NAME','obs_freq'))
        wrmean_obs = np.nanmean(wrfreq_obs,axis=1)
        pwr.plot_wr_frequency(wrmean_obs,wrnum_obs,dlen_obs,freq_plot_title,freq_plot_outname)

    if ("PLOTFREQ" in steps_list_fcst):
        # Check to be sure the forecast TIMEFREQ step has been run
        if not ("TIMEFREQ" in steps_list_fcst):
            raise Exception('Must run forecast Frequency calculation before plotting the frequencies.')
        # Get plot title, output name, compute mean, and plot
        freq_plot_title = os.environ.get('FCST_FREQ_PLOT_TITLE','Seasonal Cycle of WR Days/Week')
        freq_plot_outname = os.path.join(oplot_dir,os.environ.get('FCST_FREQ_PLOT_OUTPUT_NAME','fcst_freq'))
        wrmean_fcst = np.nanmean(wrfreq_fcst,axis=1)
        pwr.plot_wr_frequency(wrmean_fcst,wrnum_fcst,dlen_fcst,freq_plot_title,freq_plot_outname)


if __name__ == "__main__":
    main()
