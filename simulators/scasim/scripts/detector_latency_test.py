#!/usr/bin/env python
#
# Script `detector_latency_test` demonstrates how to use the
# detector latency simulation in SCASim. It can also be used to
# investigate those effects.
#
# :History:
# 
# 09 Jun 2014: Created
# 25 Jun 2014: Documentation updated and tests adjusted prior to release.
# 01 Jul 2014: Corrected amplifier gain.
# 03 Mar 2015: Simulator simplified by removing the ability to save data
#              to legacy file formats. SUB64 changed from 64x68 to 64x72.
# 27 May 2015: Replaced pyfits with astropy.io.fits
# 08 Sep 2015: SensorChipAssembly parameters are defined using a setup
#              method. Made compatible with Python 3.
# 22 Sep 2015: The detector latency test is easier now the SCA class is
#              a singleton. Added a full-frame singleton test.
# 04 Dec 2015: poisson_integrator module renamed to integrators.
#
# @author: Steven Beard (UKATC)

"""

detector_latency_test: Demonstrates SCASim detector latency simulation.

Run this script to try out the detector latency effects in SCASim.
The script may be edited to try out different simulation parameters.

There are no command line arguments.

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import numpy as np
import astropy.io.fits as pyfits

# Import MIRI data products
from miri.datamodels.miri_measured_model import MiriSlopeModel
from miri.datamodels.miri_illumination_model import \
    MiriIlluminationModel

# Import MIRI plotting module
import miri.tools.miriplot as mplt

# Import MIRI simulation utilities.
from miri.simulators.integrators import linear_regression, \
    ImperfectIntegrator
from miri.simulators.scasim.sensor_chip_assembly import SensorChipAssembly, \
    simulate_sca, simulate_sca_list

def ramps_to_slopes( array4d, grptime=1.0, diff_only=False ):
    """
    
    Convert a 4-D input array into a 3-D output array by combining
    all the group planes into a single slope plane.
    
    NOTE: This function is very inefficient!
    
    Dependency: Assumes the function linear_regression has been imported.
    
    """
    array4d = np.asarray(array4d)
    assert len(array4d.shape) == 4
    (nints, ngroups, rows, columns) = array4d.shape
    assert ngroups > 1
    array3d = np.zeros([nints, rows, columns], dtype=np.float32)
    
    if diff_only:
        # A quick and dirty method which subtracts the last ramp
        # from the first
        timediff = grptime * (ngroups - 1)
        array3d = (array4d[:,-1,:,:] - array4d[:,0,:,:]) / float(timediff)
        
    else:
        # A full slope derived from linear regression.
        timearray = grptime * np.array( list(range(0, ngroups)) )
        for intg in range(0, nints):
            for row in range(0, rows):
                for column in range(0, columns):
                    (slope, ic) = linear_regression( timearray,
                                                     array4d[intg,:,row,column])
                    array3d[intg,row,column] = slope
    return array3d

# Main program starts here.
if __name__ == "__main__":

    # There are 4 levels at which the detector latency algorithms
    # can be used and tested.
    TEST_WITH_INTEGRATOR = True
    TEST_WITH_SCA = True
    TEST_WITH_SCRIPTS = True
    TEST_FULL_FRAME = True
    
    # Define the amount of data saving and plotting
    PLOT_ILLUMINATION_DATA = False
    SAVE_ILLUMINATION_DATA = False
    
    PLOT_SLOPE_DATA = False
    SAVE_SLOPE_DATA = False
    PLOT_SEQUENCE = True

    # These parameters can be matched to JPL test data.
    NEXPOSURES = 20
    NDARKEXP = 5
    NINTS = 5
    NGROUPS = 50
    FRAME_TIME = 2.785  # Can define full-frame time for subarray-sized data.

    detectorid = 'MIRIMAGE'
    readout_mode = 'FAST'
    subarray = 'SUB64' # 'MASK1065'
    inttime = None
    det_temperature = 6.7
 
#     # Define 224 x 288 test data (MASK1065 subarray)
#     ROWS = 224
#     COLUMNS = 288
#     SAMPLE_ROW = 65
#     SAMPLE_COL = 65
#     leftdata = np.zeros([ROWS,COLUMNS])
#     leftdata[64:192,60:70] = 1.0
#     #mplt.plot_image(leftdata, title='leftdata')
#       
#     rightdata = np.zeros([ROWS,COLUMNS])
#     rightdata[64:192,188:198] = 1.0
#     #mplt.plot_image(rightdata, title='rightdata')

    # Define 64 x 72 test data (SUB64 subarray)
    ROWS = 64
    COLUMNS = 72
    SAMPLE_ROW = 16
    SAMPLE_LCOL = 16
    SAMPLE_RCOL = 48
    leftdata = np.zeros([ROWS,COLUMNS])
    leftdata[10:50,14:20] = 1.0
    if PLOT_ILLUMINATION_DATA:
        mplt.plot_image(leftdata, title='leftdata')
      
    rightdata = np.zeros([ROWS,COLUMNS])
    rightdata[10:50,44:50] = 1.0
    if PLOT_ILLUMINATION_DATA:
        mplt.plot_image(rightdata, title='rightdata')
        mplt.close()

# NOTE: FULL FRAME DATA SIMULATIONS USE A LOT OF MEMORY - MAY NEED TO LIMIT NGROUPS 
#     # Define full frame test data
#     ROWS = 1024
#     COLUMNS = 1024
#     leftdata = np.zeros([ROWS,COLUMNS])
#     r1 = 32
#     c1 = 60
#     for cc in range(c1, COLUMNS, 256):
#         for rr in range(r1, ROWS, 256):
#             leftdata[rr:rr+192,cc:cc+16] = 1.0    
#     if PLOT_ILLUMINATION_DATA:
#         mplt.plot_image(leftdata, title='leftdata')
#      
#     rightdata = np.zeros([ROWS,COLUMNS])
#     r1 = 32
#     c1 = 188
#     for cc in range(c1, COLUMNS, 256):
#         for rr in range(r1, ROWS, 256):
#             rightdata[rr:rr+192,cc:cc+16] = 1.0    
#     if PLOT_ILLUMINATION_DATA:
#         mplt.plot_image(rightdata, title='rightdata')

    # Save the illumination data to test files.
    if SAVE_ILLUMINATION_DATA:
        left_filename = "left_illumination.fits"
        right_filename = "right_illumination.fits"
        with MiriIlluminationModel(intensity=leftdata) as leftmodel:
            print( "\nTest data with left side illuminated" )
            leftmodel.set_instrument_metadata(detector=detectorid,
                                modelnam='JPL',
                                detector_temperature=det_temperature)
            leftmodel.set_exposure_metadata(readout_mode, frame_time=FRAME_TIME,
                                            nints=None, ngroups=None)
            leftmodel.set_subarray_metadata(subarray)
            print( leftmodel )
            print( "Saving to", left_filename )
            leftmodel.save(left_filename, overwrite=True)
            if PLOT_ILLUMINATION_DATA:
                leftmodel.plot()
            del leftmodel
        with MiriIlluminationModel(intensity=rightdata) as rightmodel:
            print( "\nTest data with right side illuminated" )
            rightmodel.set_instrument_metadata(detector=detectorid,
                                modelnam='JPL',
                                detector_temperature=det_temperature)
            rightmodel.set_exposure_metadata(readout_mode, frame_time=FRAME_TIME,
                                             nints=None, ngroups=None)
            rightmodel.set_subarray_metadata(subarray)
            print( rightmodel )
            print( "Saving to", right_filename )
            rightmodel.save(right_filename, overwrite=True)
            if PLOT_ILLUMINATION_DATA:
                rightmodel.plot()
            del rightmodel

    # These multipliers are used to define the various levels of
    # illumination described in the JPL 3 tests.
    dark_dn = 0.0
    very_faint_dn = 0.5
    faint_dn = 2.0
    medium_faint_dn = 5.0
    medium_dn = 20.0
    med_bright_dn = 80.0
    bright_dn = 400.0
    very_bright_dn = 1000.0

    ampgain = 0.1818      # DN/s = electrons/s * ampgain

    # Define the sequence of observations to be simulated.
    # Instead of varying the brightness of the whole array, this sequence
    # illuminates the left and right sides of the array with a different
    # sequence.
    scale_sequence = [med_bright_dn, medium_dn,  med_bright_dn, medium_dn]
    name_sequence  = ['DARK-MBRI',   'MED-DARK', 'MBRI-DARK',  'DARK-MED']
    data_sequence  = ['R',           'L',        'L',          'R']

    #persistence = [1.0e-8, 0.03, 0.0]
    persistence = None
    #linearity = [1.0, 0.0]
    linearity = None
    zpslow = [50000.0, 0.0084]
    zpfast = [[0.0, -2.917],
              [0.0, -2.292],
              [0.0, -2.396],
              [0.0, -2.408]]
    slow_gain =  1.67e-9   # Scale factor for slow latent
    slow_decay = 136000.0  # Timescale of slow latent decay (s)
    fast_gain = 0.002      # Scale factor for fast latent
    fast_decay = 300.0     # Timescale of fast latent decay (s)

    #
    # (1) TEST USING INTEGRATOR OBJECT DIRECTLY (NO AMPLIFIERS OR REF PIXELS)
    #
    if TEST_WITH_INTEGRATOR:
        print( "\n(1) Testing detector latency effects with integrator object" )
        integrator = ImperfectIntegrator(ROWS, COLUMNS, bucket_size=250000,
                                         verbose=4)
        integrator.simulate_poisson_noise = False
        integrator.set_persistence(persistence)
        integrator.set_linearity(linearity)
        integrator.set_zeropoint(zpslow, zpfast)
        integrator.set_latency([slow_gain,slow_decay], [fast_gain,fast_decay])
        latency_string = "%s; %s" % (integrator._slow_latency_str(),
                                     integrator._fast_latency_str())
        print( integrator )

        time = []
        ctime = []
        lsignal = []
        rsignal = []
        lslope = []
        rslope = []
        lrelslope = []
        rrelslope = []
        elapsed = 0.0
 
        sno = 0
        for sflux in scale_sequence:
            snono = sno + 1
            sname = name_sequence[sno]
            if data_sequence[sno] == 'L':
                fluxdata = leftdata * sflux / ampgain
            else:
                fluxdata = rightdata * sflux / ampgain
            for exp in range(0,NEXPOSURES):
                fullarray = np.zeros([NINTS,NGROUPS,ROWS,COLUMNS])
                expno = exp+1
                title = "\n%s source (%.2f DN/s), exposure %d" % (sname, sflux, expno)
                print( title )
                for integ in range(0,NINTS):
                    print( "Integration", integ+1 )
                    if integ == 0:
                        # There are more resets before the first integration
                        # of a new exposure.
                        integrator.reset(nresets=3, new_exposure=True)
                    else:
                        integrator.reset(new_exposure=False)
                    ctime.append(elapsed)
                    for group in range(0,NGROUPS):
                        integrator.integrate(fluxdata, FRAME_TIME)
                        fullarray[integ,group,:,:] = integrator.readout()
                        elapsed += FRAME_TIME
                        time.append(elapsed)
                        lsignal.append(fullarray[integ,group,SAMPLE_ROW,
                                                 SAMPLE_LCOL])
                        rsignal.append(fullarray[integ,group,SAMPLE_ROW,
                                                 SAMPLE_RCOL])
                        
                slope_data = ramps_to_slopes( fullarray, grptime=FRAME_TIME,
                                              diff_only=False )
                slope_data *= ampgain # Convert back to DN/s
                print( "Slope data range: %g to %g" % \
                    (slope_data.min(), slope_data.max()) )
                if exp == 0:
                    firstslope = slope_data
                for integ in range(0,NINTS):
                    lslope.append(slope_data[integ, SAMPLE_ROW,SAMPLE_LCOL])
                    rslope.append(slope_data[integ, SAMPLE_ROW,SAMPLE_RCOL])
                    sdiff = slope_data[integ, SAMPLE_ROW,SAMPLE_LCOL] - \
                                firstslope[0, SAMPLE_ROW,SAMPLE_LCOL]
                    lrelslope.append(sdiff)
                    sdiff = slope_data[integ, SAMPLE_ROW,SAMPLE_RCOL] - \
                                firstslope[0, SAMPLE_ROW,SAMPLE_RCOL]
                    rrelslope.append(sdiff)
                # Only save and display a selection of the exposures.
                if (expno < 3) or ((expno % 10) == 0):
                    if PLOT_SLOPE_DATA:
                        mplt.plot_image(np.abs(slope_data), datascale='linear',
                                    zlabel='Integration', title=title)
                    if SAVE_SLOPE_DATA:
                        filename = "1_detector_latency%.02d_%.02d.fits" % (snono,expno)
                        with MiriSlopeModel( data=slope_data ) as slope_model:
                            #slope_model.update( sca.exposure_data ) # Update metadata
                            slope_model.save(filename, overwrite=True)
                            del slope_model
                del slope_data
                del fullarray
            sno += 1
            del firstslope
            del fluxdata
            if PLOT_SLOPE_DATA:
                mplt.close()

        # Display ramp and slope plots for the whole observation sequence.
        print( "There are", len(time), "ramps" )
        if PLOT_SEQUENCE:
            ramps_title = "Ramps for sequence %s" % str(name_sequence)
            ramps_title += "\n" + latency_string
            mplt.plot_xycolumn(time, (lsignal,rsignal), figsize=(20,8),
                           xlabel='Time (seconds)',
                           ylabels=('Left reading DN', 'Right reading (DN)'),
                           title=ramps_title)
        print( "There are", len(ctime), "slopes" )
        if PLOT_SEQUENCE:
            slopes_title = "Slopes for sequence %s" % str(name_sequence)
            slopes_title += "\n" + latency_string
            mplt.plot_xycolumn(ctime, (lslope,rslope), figsize=(20,8),
                           xlabel='Time (seconds)',
                           ylabels=('Left slope (DN/s)', 'Right slope (DN/s)'),
                           title=slopes_title)
            rel_title = "Slope relative to first %s" % str(name_sequence)
            rel_title += "\n" + latency_string
            mplt.plot_xycolumn(ctime, (lrelslope,rrelslope), figsize=(20,8),
                           xlabel='Time (seconds)',
                           ylabels=('Left subtracted slope (DN/s)',
                                    'Right subtracted slope (DN/s)'),
                           title=rel_title)
            mplt.close()

    # These parameters control the SCASim simulation
    cosmic_ray_mode = 'NONE'
    qe_adjust = False
    simulate_poisson_noise=False
    simulate_read_noise=False
    simulate_ref_pixels=False
    simulate_bad_pixels=False
    simulate_dark_current=False
    simulate_amp_effects=False

    #
    # (2) TEST USING SENSOR CHIP ASSEMBLY OBJECT DIRECTLY
    #
    if TEST_WITH_SCA:
        print( "\n(2) Testing detector latency effects with SCA object" )
     
        # Create a sensor chip assembly object with the required readout mode.
        sca = SensorChipAssembly()
        sca.setup(detectorid, readout_mode=readout_mode,
                  subarray=subarray, inttime=inttime, ngroups=NGROUPS,
                  nints=NINTS, temperature=det_temperature,
                  cosmic_ray_mode=cosmic_ray_mode,
                  qe_adjust=qe_adjust,
                  simulate_poisson_noise=simulate_poisson_noise,
                  simulate_read_noise=simulate_read_noise,
                  simulate_ref_pixels=simulate_ref_pixels,
                  simulate_bad_pixels=simulate_bad_pixels,
                  simulate_dark_current=simulate_dark_current,
                  simulate_amp_effects=simulate_amp_effects,
                  makeplot=False, verbose=3)
     
        sno = 0
        for sflux in scale_sequence:
            snono = sno + 1
            sname = name_sequence[sno]
            if data_sequence[sno] == 'L':
                fluxdata = leftdata * sflux / ampgain
            else:
                fluxdata = rightdata * sflux / ampgain
            illumination_map = MiriIlluminationModel(intensity=fluxdata)
            sca.set_illumination(illumination_map)
            # Override the properties defined in detector_properties.py
            # NOTE: The sca.detector.pixels object is the same as the
            # integrator object in test (1).
            sca.detector.pixels.set_persistence(None)
            sca.detector.pixels.set_linearity(None)
            sca.detector.pixels.set_zeropoint(zpslow, zpfast)
            sca.detector.pixels.set_latency([slow_gain,slow_decay],
                                            [fast_gain,fast_decay])
            #print( sca )
            for exp in range(0,NEXPOSURES):
                expno = exp+1
                title = "\n%s source (%.2f DN/s), exposure %d" % (sname, sflux, expno)
                print( title )
                simdata = sca.exposure(frame_time=FRAME_TIME)
     
                # Only save and display a selection of the exposures.
                if (expno < 3) or ((expno % 10) == 0):
                    #sca.plot()
                    filename = "2_detector_latency%.02d_%.02d.fits" % (snono,expno)
                    #sca.write_data(filename, overwrite=True)
#                     firstint = sca.exposure_data.get_integration(0)
#                     leftramp = firstint[:, SAMPLE_ROW, SAMPLE_LCOL]
#                     print( "leftramp=", leftramp )
#                     sca.plot_ramp(SAMPLE_ROW, SAMPLE_LCOL,
#                                   description=title + ' L',
#                                   frame_time=FRAME_TIME, show_ints=True)
#                     rightramp = firstint[:, SAMPLE_ROW, SAMPLE_RCOL]
#                     print( "rightramp=", rightramp )
#                     sca.plot_ramp(SAMPLE_ROW, SAMPLE_RCOL,
#                                   description=title + 'R',
#                                   frame_time=FRAME_TIME, show_ints=True)
                    slope_data = sca.exposure_data.slope_data(diff_only=False)
                    if not simulate_amp_effects:
                        slope_data *= ampgain # Convert back to DN/s
                    print( "Slope data range: %g to %g" % \
                        (slope_data.min(), slope_data.max()) )
                    if PLOT_SLOPE_DATA:
                        mplt.plot_image(slope_data, datascale='linear',
                            zlabel='Integration', withbar=False,
                            title='Slope data for %s source exposure %d' % \
                            (sname, expno))
                    if SAVE_SLOPE_DATA:
                        with MiriSlopeModel( data=slope_data ) as slope_model:
                            #slope_model.update( sca.exposure_data ) # Update metadata
                            slope_model.save(filename, overwrite=True)
                            del slope_model
                    del slope_data
                del simdata
            sno += 1
            if PLOT_SLOPE_DATA:
                mplt.close()
        del sca
    
    # (3) TEST USING HIGH LEVEL SCASIM SCRIPTS
    if TEST_WITH_SCRIPTS:
        print( "\n(3) Testing detector latency effects with SCASim scripts" )
        
        # This test will only work if there are some illumination
        # data files available.
        if SAVE_ILLUMINATION_DATA:
            
            # Set up a list of input and output file names
            inputfiles = [left_filename, right_filename,
                          left_filename, right_filename]
            outputfiles = ['LatentTestOutput1.fits', 'LatentTestOutput2.fits',
                           'LatentTestOutput3.fits', 'LatentTestOutput4.fits']
            
            simulate_sca_list(inputfiles, outputfiles,
                    detectorid, readout_mode=readout_mode,
                    subarray=subarray, inttime=inttime, ngroups=NGROUPS,
                    nints=NINTS, temperature=det_temperature,
                    cosmic_ray_mode=cosmic_ray_mode,
                    qe_adjust=qe_adjust,
                    simulate_poisson_noise=simulate_poisson_noise,
                    simulate_read_noise=simulate_read_noise,
                    simulate_ref_pixels=simulate_ref_pixels,
                    simulate_bad_pixels=simulate_bad_pixels,
                    simulate_dark_current=simulate_dark_current,
                    simulate_amp_effects=simulate_amp_effects,
                    overwrite=True, makeplot=False, verbose=2)
            
        else:
            print( "Test skipped - no illumination data files available." )

    # (4) FULL FRAME DECAY TEST
    if TEST_FULL_FRAME:
        print( "\n(4) Full frame decay test with Singleton" )
        
        # This test will only work if there are some data files available.
        flashfile = "TestOutputBOX.fits"
        decayfile = "TestOutputIMAGE.fits"
        outputstub = "LatentDecayOutput"
        if os.path.isfile(flashfile) and os.path.isfile(decayfile):
            
            # Simulate flashing the detectors with a bright box
            ii = 1
            outputfile = "%s%.2d.fits" % (outputstub, ii)
            simulate_sca(flashfile, outputfile,
                    detectorid, readout_mode=readout_mode,
                    subarray='FULL', inttime=inttime, ngroups=10,
                    nints=NINTS, temperature=det_temperature,
                    cosmic_ray_mode=cosmic_ray_mode,
                    qe_adjust=qe_adjust,
                    simulate_poisson_noise=simulate_poisson_noise,
                    simulate_read_noise=simulate_read_noise,
                    simulate_ref_pixels=simulate_ref_pixels,
                    simulate_bad_pixels=simulate_bad_pixels,
                    simulate_dark_current=simulate_dark_current,
                    simulate_amp_effects=simulate_amp_effects,
                    overwrite=True, makeplot=False, verbose=2)
            
            # Simulate observing test data.
            for ii in (2,3,4):
                outputfile = "%s%.2d.fits" % (outputstub, ii)
                simulate_sca(decayfile, outputfile,
                    detectorid, readout_mode=readout_mode,
                    subarray='FULL', inttime=inttime, ngroups=10,
                    nints=NINTS, temperature=det_temperature,
                    cosmic_ray_mode=cosmic_ray_mode,
                    qe_adjust=qe_adjust,
                    simulate_poisson_noise=simulate_poisson_noise,
                    simulate_read_noise=simulate_read_noise,
                    simulate_ref_pixels=simulate_ref_pixels,
                    simulate_bad_pixels=simulate_bad_pixels,
                    simulate_dark_current=simulate_dark_current,
                    simulate_amp_effects=simulate_amp_effects,
                    overwrite=True, makeplot=False, verbose=2)        
        else:
            print( "Test skipped - no illumination data files available." )

    print( "Tests finished." )