#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module spec_tools - Contains general purpose functions for analysing MIRI
spectroscopic data, including:

condense - function to rebin, smooth, filter (spectrally) 1 dimensional arrays

convGauss - smoothing of a spectrum array by convolving with a Gaussian

lrs2D_spextract -Extracts a 2-D spectrum from 2-D data by applying a simple
    rectangular aperture

lrs_extract_spec - Extracts spectra from the data by summing along the
    spatial-x axis, considering NaNs.

lrs_extract_spec_with_fit - As above, but the row extension to sum over is
    determined by a Gaussian fit +- 3*FWHM

get_psf_fit - Fits a pinhole (point source) continuum spectrum by a Gaussian
    for each row (x axis).

optimalSpecExtraction - Extracts spectra from the data by summing along the
    spatial-x axis considering NaNs and error weights.

subtractBackground - Subtracts the background from a spectrum, considers and
    avoids the latents (transients) by taken the median of the background, 
     and propagates the errors.

lrs_round - round like in python2 built-in rounding way

interpolWaveOnRows - interpolates the wavelengths and positions arrays from the wavelength calibration
    file to the detector rows

interpolLin - linear interpolation of a spectrum to new wavelengths

interpolSpline - B-spline interpolation of a spectrum to new wavelengths

calcStraightLine - Calculates the slope and the axis intercept of a straight
    line using points in an x, y grid


Longer description, references, class hierarchy...

:History:

20 Sep 2012: Created
13 Nov 2012: Major restructuring of the package folder.
             Import statements updated.
15 Mar 2013: condense and convGauss created
12 Jul 2013: Original library of spectral flattening tools.
20 Feb 2014: Old data product removed.
02 Dec 2014: Replaced DQ max function to prevent failure
             of unit test.
09 Dec 2014: Suppress unwanted output.
15 Jun 2017: Corrected issue with non-integer numpy indices (which fail
             under numpy 1.12.1).
08 Nov 2017: Moved from LRS pipeline package to general purpose spec_tools
             package.
17 Oct 2019: Merged with the spectral flattening tools from MiriTeam

@author:  Juergen Schreiber, Steven Beard (UKATC)

"""

# Import standard python packages/modules
# e.g. import os

# Import packages/modules from the MIRI development environment
import numpy as np
# "numba" import disabled due to segfault issue on import when Numba has been
# installed via pip into a Conda environment, see MIRI-749.
# import numba as nb
import scipy.signal as sig
from scipy import interpolate

# Import MIRI packages/modules with full namespace
from miri.datamodels.miri_measured_model import MiriMeasuredModel
#from miri.spectroscopy import Spectrum1D
from miri.tools.fitting import nonLinFit, gaussian
from matplotlib import pyplot
import miri.tools.miriplot as mplt
import math

def condense(spec, chunk, method = "Median", running = True):
    """
    
    smoothing/binning of an spectrum array
    
    :Parameters:
    
    spec: array 
        spectrum
    
    chunk: int 
        number of elements to average
    
    method: string, optional (default = "Median") 
        method of smoothing /rebinning: Mean or Median or Min
    
    running: boolean, optional (default = True)
        true if its a running averaging (smoothing) (default),
        false if it is rebinning 
   
    :Returns:
    
    rebinned/smoothed spectrum array
    
    """        
    if not isinstance(chunk, int):
        print("in condense: chunk parameter is not integer, will be rounded to integer:")
        chunk = int(np.rint(chunk))
        print(str(chunk))
    leng = len(spec)
    if leng < chunk:
        print("no of elements in spec is smaller than chunk value!")
        raise ValueError
    no_chunks = np.int32(leng/chunk)
    newspec = np.copy(spec)
    if running:
        newspec = np.copy(spec)
        for i in range(leng):
            if method == "Median":   
                if i + chunk < leng:
                    newspec[i] = np.median(spec[i:i+chunk])
                else:
                    newspec[i] = np.median(spec[i:])
            elif method == "Mean":
                if i + chunk < leng:
                    newspec[i] = np.mean(spec[i:i+chunk])
                else:
                    newspec[i] = np.mean(spec[i:])
            elif method == "Min":
                if i + chunk < leng:
                    newspec[i] = np.min(spec[i:i+chunk])
                else:
                    newspec[i] = np.min(spec[i:])
            else:
                raise ValueError("method parameter has wrong value") 
    else:
        newspec = np.zeros(no_chunks)
        for i in range(no_chunks):
            if method == "Median":
                if (i+1) * chunk < leng:
                    newspec[i] = np.median(spec[i*chunk:(i+1)*chunk])
                else:
                    newspec[i] = np.median(spec[i*chunk:])
                    break
            elif method == "Mean":
                if (i+1) * chunk < leng:
                    newspec[i] = np.mean(spec[i*chunk:(i+1)*chunk])
                else:
                    newspec[i] = np.mean(spec[i*chunk:])
                    break 
            elif method == "Min":
                if i + chunk < leng:
                    newspec[i] = np.min(spec[i*chunk:(i+1)*chunk])
                else:
                    newspec[i] = np.min(spec[i*chunk:])
                    break            
            else:
                raise ValueError               
            
    return newspec

def convGauss(x, n_points, sigma, mode = 'valid'):
    """
    
    smoothing of an spectrum array by convolving with a Gaussian
    
    :Parameters:
    
    x: array 
        spectrum
    
    n_points: int
        number of elements of Gaussian kernel
    
    sigma: float
        The stdev, sigma of Gaussian function
    
    mode: str, optional (default = "valid")
        {'full', 'valid', 'same'}
            * 'full':
                This returns the convolution at each point of overlap, with an output
                shape of (N+M-1,). At the end-points of the convolution, the signals
                do not overlap completely, and boundary effects may be seen.
            * 'same':
                Mode same returns output of length max(M, N). Boundary effects are
                still visible.
            * 'valid': (default)
                Mode valid returns output of length max(M, N) - min(M, N) + 1.
                The convolution product is only given for points where the signals
                overlap completely. Values outside the signal boundary have no effect.
    
    :Returns:
    
    smoothed spectrum array
    
    """         
    gauss = sig.gaussian(n_points, sigma)
    return np.convolve(x, gauss/np.sum(gauss), mode)


def lrs2D_spextract(data, xmin = 0, xmax = 1024, ymin =0, ymax = 1032, copy = True):
    """

    Extracts a 2-D spectrum from the 2-D data product of the imager by applying a simple
    rectangular aperture.

    :Parameters:

    data: MiriMeasuredModel (or similar 2-D DataModel)
        The 2-D data product from which the 1-D spectrum will be extracted.
        
    xmin, xmax: int, optional (default=0, 1024)
        aperture x coordinates,
        The aperture from which the spectrum is to be extracted

    ymin, ymax: int, optional (default=0, 1032) 
        aperture y coordinates,
        The aperture from which the spectrum is to be extracted

    :Returns:
    
    array: LRS Measured data product 
        The extracted 2D area of the imager product

    :Examples:

    imager_product = MiriMeasuredModel( filename )
    
    lrs_product = lrs_2Dspectract( imager_product )

    """
    # Make sure the given data is a 2-D JWST data product.
    if isinstance(data, MiriMeasuredModel):
        if np.ndim(data.data) != 2:
            strg = "lrs2D_spextract - 2 dimensional data expected "
            strg += "(%d given)." % np.ndim(data.data)
            raise TypeError(strg)
    else:
        raise TypeError("lrs2D_spextract - " + \
                        "Given data product is not a MiriMeasuredModel")
    if np.ndim(data.data)==0:
        raise TypeError("lrs2D_spextract - " + \
                        "Given data product does not contain any data!")


    if copy:
        data = data.copy()
        
    if np.ndim(data.data) > 1:
        sciarray = data.data[xmin:xmax, ymin:ymax]
        data.data=sciarray
        #data.replace_dataarray('SCI',sciarray , dtype = np.float, copy = True)
            
    if np.ndim(data.err) > 1:
        errarray = data.err[xmin:xmax, ymin:ymax]
        data.err = errarray
        #data.replace_dataarray('ERR',errarray , dtype = np.float, copy = True)       
    
    if np.ndim(data.dq) > 1:
        dqarray = data.dq[xmin:xmax, ymin:ymax]
        data.dq = dqarray
        #data.replace_dataarray('DQ',dqarray , dtype = np.uint8, copy = True)  
    #data.update_metadata('SCI')
           
    return data


def lrs_extract_spec(data):
    """
    
    Extracts spectra from the data by summing along the spatial-x axis considering NaNs.

    :Parameters:

    data: MiriMeasuredModel (or similar 2-D DataModel)
        The 2-D data product from which spectra will be extracted.
        
    :Returns:
    
    Spectrum1d: dataproduct containing
        * spec: 1d array 
            The extracted spectrum.
        * err: 1d array
            The propagated error 
    
    """
   
    # Make sure the given data is a 2-D JWST data product.
    if isinstance(data, MiriMeasuredModel):
        if np.ndim(data.data) != 2:
            strg = "lrs_extract_spec - 2 dimensional data expected "
            strg += "(%d given)." % np.ndim(data.data)
            raise TypeError(strg)
    else:
        raise TypeError("lrs_extract_spec - " + \
                        "Given data product is not a MiriMeasuredModel")
    if np.ndim(data.data) == 0:
        raise TypeError("lrs_extract_spec - " + \
                        "Given data product does not contain any data!")
    
    da = data.data
    spec = np.nansum(da, 1)
    data.data = spec
    if np.ndim(data.err) > 0:
        er = data.err
        var = np.nansum(np.square(er),1)
        err = np.sqrt(var)
        data.err = err
    else:
        err = None
    if np.ndim(data.dq) > 1:
        dq = data.dq
        # Crunch the DQ array into 1-D by taking the maximum
#         qual = np.max(dq, 1)
        qual = dq.max(0)
        data.dq = qual
    else:
        qual = None
    
    
    
    #spectrum = Spectrum1D(spec, error = err, quality = qual )
    return data


_sigma_clip = 1.
def lrs_extract_spec_with_fit(data, minX = 0, leftY = 0, verbose = False,
                              makeplot = False):
    """
    
    Extracts spectra from the data by summing along the spatial-x axis considering NaNs.
    The row extension to sum over is determined by a Gaussian fit +- 3*FWHM

    :Parameters:

    data: MiriMeasuredModel (or similar 2-D DataModel)
        The 2-D data product from which spectra will be extracted.
        
    :Returns:
    
    Spectrum1d: dataproduct with
        * spec: 1d array 
            The extracted spectrum
        * err: 1d array
            The propagated error 

    :Raises:
    
    TypeError
        if wrong input type
    
    """
   
    # Make sure the given data is a 2-D JWST data product.
    if isinstance(data, MiriMeasuredModel):
        if np.ndim(data.data) != 2:
            strg = "lrs_extract_spec - 2 dimensional data expected "
            strg += "(%d given)." % np.ndim(data.data)
            raise TypeError(strg)
    else:
        raise TypeError("lrs_extract_spec - " + \
                        "Given data product is not a MiriMeasuredModel")
    if np.ndim(data.err) != 2:
        raise TypeError("lrs_extract_spec - " + \
                        "Given data product does not contain any data!")
    
    sig = data.data
    sigErr = data.err
    col = []
    xpixel = []
    fwhm = []
    for i in range(sig.shape[0]):
        si = sig[i,:]
        simax = np.max(si)
        stdevSi = np.std(si)
        if simax <= 0 or stdevSi == 0:
            continue
        if simax < 0.03 and simax/stdevSi < 1.:
            continue
        maxind = np.where(si == simax)[0][0]
        p_guess = (0, simax, maxind,2)
        try:
            fit, p, perr = nonLinFit(gaussian, np.arange(sig.shape[1]), si, y_sigma = sigErr[i,:], p_guess = p_guess, plotting = False, verbose=verbose)
            col.append(p[2] + leftY)
            xpixel.append(i+minX)
            fwhm.append(p[3])
        except:
            print("fit ",i, " did not work")
        
    col = np.asarray(col)
    xpixel = np.asarray(xpixel)
    fwhm = np.asarray(fwhm)
    if makeplot:
        mplt.plot_xy(xpixel,col, title="fit result position of peak")
    stdevCol = np.std(col)
    meanCol = np.mean(col)
    if verbose:
        print("mean column ", meanCol)
        print("stddev column ", stdevCol)
    colInd = np.where(np.logical_and(col < meanCol + _sigma_clip*stdevCol, col > meanCol - _sigma_clip*stdevCol))[0]
        
    newXpixel = xpixel[colInd]
    newCol = col[colInd]
    newFwhm = fwhm[colInd]
    coeff = np.polyfit(newXpixel, newCol, 1)
    line = coeff[0] * newXpixel + coeff[1]
    fwhmCoeff = np.polyfit(newXpixel, newFwhm, 1)
    lineFwhm = fwhmCoeff[0] * newXpixel + fwhmCoeff[1]
    if verbose:
        print("no of useful pixels ", len(newXpixel))
    if makeplot:            
        fig = pyplot.figure(facecolor="0.98")
        fit = fig.add_subplot(211)
        fit.plot(newXpixel,newCol,'ro', np.sort(newXpixel), line)
        fit.set_title("fitted determined columns")
        fw = fig.add_subplot(212)
        fw.plot(newXpixel, lineFwhm)
        fw.set_title("FWHM")
        pyplot.show()
       
                
    #extract the spectra along col +- 3 * FHWM
    spec = np.ndarray(sig.shape[0])
    err = np.ndarray(sig.shape[0]) 
    for i in range(sig.shape[0]):
        column = coeff[0] * (i+minX) + coeff[1] - leftY
        fw = fwhmCoeff[0] * (i+minX) + fwhmCoeff[1]  
        si = sig[i, int(column - 3. * fw) : int(column + 3. * fw)]
        er = sigErr[i, int(column - 3. * fw) : int(column + 3. * fw)] 
        spec[i] = np.nansum(si)
        err[i] = np.sqrt(np.nansum(np.square(er)))
    
    if np.ndim(data.dq) > 1:
        dq = data.dq
        # Crunch the DQ array into 1-D by taking the maximum
#         qual = np.max(dq, 1)
        qual = dq.max(0)
        data.dq = qual
    else:
        qual = None
    
    data.data = spec
    data.err = err
    #spectrum = Spectrum1D(spec, error = err, quality = qual )
    return data, newXpixel, line


def get_psf_fit(data, minX = 0, leftY = 0, verbose = False, makeplot = False):
    """
    
    fits a pinhole (point source) continuum spectrum by a Gaussian for each row (x axis).
    
    :Parameters:

    data: MiriMeasuredModel (or similar 2-D DataModel)
        The 2-D data product from which spectra will be extracted.
        
    :Returns:
    
    * x values per row
    * fit function (vs x) per row
    * all Gaussian fit parameters per row
    * all errors of Gaussian fit parameters per row
    
    """
   
    # Make sure the given data is a 2-D JWST data product.
    if isinstance(data, MiriMeasuredModel):
        if np.ndim(data.data) != 2:
            strg = "lrs_extract_spec - 2 dimensional data expected "
            strg += "(%d given)." % np.ndim(data.data)
            raise TypeError(strg)
    else:
        raise TypeError("lrs_extract_spec - " + \
                        "Given data product is not a MiriMeasuredModel")
    if np.ndim(data.err) != 2:
        raise TypeError("lrs_extract_spec - " + \
                        "Given data product does not contain any data!")
    
    sig = data.data
    sigErr = data.err
    all_x = []
    all_p = []
    all_fit = []
    all_perr = []
    for i in range(sig.shape[0]):
        si = sig[i,:]
        simax = np.max(si)
        stdevSi = np.std(si)
        if simax <= 0 or stdevSi == 0:
            if verbose:
                print("signal is below or equal zero in row ", i)
            continue
        if simax/stdevSi < 1. or simax < np.mean(sigErr):
            if verbose:
                print("signal too low in row ", i)
            continue
        maxind = np.where(si == simax)[0][0]
        p_guess = (0, simax, maxind,2)
        try:
            fit, p, perr = nonLinFit(gaussian, np.arange(sig.shape[1]), si, y_sigma = sigErr[i,:], p_guess = p_guess, plotting = makeplot, verbose=verbose)
            p[2] = p[2]+leftY
            all_p.append(p)
            all_fit.append(fit)
            all_perr.append(perr)
            all_x.append(i + minX)
        except:
            print("fit ",i, " did not work")
                
    return all_x, all_fit, all_p, all_perr


def optimalSpecExtraction(data):
    """
    
    Extracts spectra from the data by summing along the spatial-x axis considering NaNs and error weights.

    :Parameters:

    data: MiriMeasuredModel (or similar 2-D DataModel)
        The 2-D data product from which spectra will be extracted.
        
    :Returns:
    
    Spectrum1d: dataproduct with
        spec: 1d array 
            The extracted spectrum.
        err: 1d array
            The propagated error 
    
    """
   
    # Make sure the given data is a 2-D JWST data product.
    if isinstance(data, MiriMeasuredModel):
        if np.ndim(data.data) != 2:
            strg = "optimalSpecExtraction - 2 dimensional data expected "
            strg += "(%d given)." % np.ndim(data.data)
            raise TypeError(strg)
    else:
        raise TypeError("optimalSpecExtraction - " + \
                        "Given data product is not a MiriMeasuredModel")

    if np.ndim(data.err) != 2:
        raise TypeError("optimalSpecExtraction - " + \
                        "Given data product does not contain errors, optimal extraction is not possible!")    
    da = data.data
    er = data.err
    var = np.square(er)
    weight = 1./var
    sum_weight = np.nansum(weight, 1)
    
    lenX = len(sum_weight)
    spec = np.zeros(lenX)
    for i in range(lenX):
        spec[i] = np.nansum(da[i, :] * weight[i, :]/sum_weight[i])
    
    err = np.sqrt(1./sum_weight)

    data.data = spec
    data.err = err
    if data.has_dataarray('DQ'):
        dq = data.dq
        # Crunch the DQ array into 1-D by taking the maximum
#         qual = np.max(dq, 1)
        qual = dq.max(0)
        data.dq = qual
    else:
        qual = None
    
    #spectrum = Spectrum1D(spec, error = err, quality = qual )
    return data


def subtractBackground(on, background):
    """
    
     subtract the background from the on spectrum and return a MiriMeasuredModel
     consider and avoid the latents (transients) by taken the median of the background, 
     and propagates the errors
     
    :Parameters:

    on: MiriMeasuredModel (or similar 2-D DataModel)
        The 2-D data product of the on measurement, LRS area already cut out

    background: MiriMeasuredModel (or similar 2-D DataModel)
        The 2-D data product of the background measurement, LRS area already cut out
    
 
    :Returns:
    
    MiriMeasuredModel 
        with the same dimensions
     
    """
    # Make sure the given data is a 2-D JWST data product.
    if isinstance(on, MiriMeasuredModel):
        if np.ndim(on.data) != 2:
            strg = "subtractBackground - 2 dimensional data expected "
            strg += "(%d given)." % np.ndim(on.data)
            raise TypeError(strg)
    else:
        raise TypeError("subtractBackground - " + \
                        "Given data product is not a MiriMeasuredModel")
    if np.ndim(on.data) == 0:
        raise TypeError("subtractBackground - " + \
                        "Given data product does not contain any data!")

    if isinstance(background, MiriMeasuredModel):
        if np.ndim(background.data) != 2:
            strg = "subtractBackground - 2 dimensional data expected "
            strg += "(%d given)." % np.ndim(background.data)
            raise TypeError(strg)
    else:
        raise TypeError("subtractBackground - " + \
                        "Given data product is not a MiriMeasuredModel")
    if np.ndim(background.data) == 0:
        raise TypeError("subtractBackground - " + \
                        "Given data product does not contain any data!")


    onData = on.data
    backgroundData = background.data
    backgroundErr = background.err
    if onData.shape[0] != backgroundData.shape[0] or onData.shape[1] != backgroundData.shape[1]:
        raise TypeError("dimensions of the data arrays of on and background are different")

    newBack = np.median(backgroundData,1)
    newBackErr = np.median(backgroundErr,1)
    # copy 1d array to 2d
    lenX = backgroundData.shape[1]
    for i in range(lenX):
        backgroundData[:,i] = newBack
        backgroundErr[:,i] = newBackErr

    background.data = backgroundData
    background.err  = backgroundErr    
           
    return on - background   
   

def lrs_round(floatNumber):
    """
    rounding floats to integers in a way that python 2 did
    e.g. lrs_round(4.5) results in 5
    
    :Parameters:
    
    floatNumber: float to be rounded to integer
    
    :Returns:
    
    rounded number  
    """
    return floatNumber/abs(floatNumber) * math.floor(abs(floatNumber)+0.5)

def interpolWaveOnRows(pos, wave):
    """
    
    interpolates the wavelengths and positions arrays from the wavelength calibration
    file to the detector rows
    
    :Parameters:
    
    pos: array 
        float positions on the array
    
    wave: array 
        corresponding wavelengths
   
    :Returns:
    
    rows and corresponding interpolated wavelengths
   
    """        
    rows = np.arange(0, lrs_round(np.nanmax(pos)), 1.)
    inter = interpolate.interp1d(pos, wave, bounds_error = False)
    new_wave= inter(rows)
    #a,b,c = interpolate.splrep(pos, wave, k=3)
    #new_wave = interpolate.splev(rows, (a,b,c))
    
    return rows, new_wave

def interpolLin(wave, spec, new_wave):
    """
    
    linear interpolation of a spectrum to new wavelengths
    
    :Parameters:
    
    wave: array 
        original wavelengths
    
    spec: array 
        spectrum
    
    new_wave: array 
        corresponding new wavelengths
   
    :Returns:
    
    spectrum: array 
        interpolated to new wavelengths
    
    """        
    inter = interpolate.interp1d(wave, spec, bounds_error = False)
    return inter(new_wave)

# Numba JIT disabled due to issue with Numba, see MIRI-749.
# @nb.jit(nopython=True, parallel=True)
def interpol_lin(wave, spec, new_wave, extra = True):
    """
    
    linear interpolation of a spectrum to new wavelengths without using scipy interpolation
    
    :Parameters:
    
    wave: array 
        original wavelengths
    
    spec: array 
        spectrum
    
    new_wave: array 
        corresponding new wavelengths

    extra: boolean
        if True, extrapolates with edge values
        if False, sets the rest to NaN
   
    :Returns:
    
    spectrum: array 
        interpolated to new wavelengths
    
    """        
    if len(wave) != len(spec):
        raise Exception("dimensions of wavelength and spectrum array are different")
    new_dim = len(new_wave)
    ret_spec = np.zeros(new_dim)
    for i in range(new_dim):
        nw = new_wave[i]
        if extra:
            if nw <= wave[0]:
                ret_spec[i]= spec[0]
            if nw >= wave[-1]:
                ret_spec[i] = spec[-1]
        else:
            if (nw < wave[0] or nw > wave[-1]):
                ret_spec[i] = np.NaN
                continue
        wave_diff_abs = np.abs(wave - nw)
        min_wave_diff = np.min(wave_diff_abs) 
        ind = np.where(wave_diff_abs == min_wave_diff)[0][0]
        wave_diff = wave[ind] - nw
        
        if wave_diff == 0:
            ret_spec[i] = spec[ind]
            continue
        elif wave_diff < 0:
            pre_ind = ind-1
            post_ind = ind
        else: 
            pre_ind = ind
            post_ind = ind+1  
        
        weight = (nw - wave[pre_ind])/(wave[post_ind]-wave[pre_ind])
        ret_spec[i] = spec[pre_ind] + weight* (spec[post_ind]-spec[pre_ind])
    
    return ret_spec

def interpolSpline(wave, spec, new_wave):
    """
    
    B-spline interpolation of a spectrum to new wavelengths
    
    :Parameters:
    
    wave: array 
        original wavelengths
    
    spec: array 
        spectrum
    
    new_wave: array 
        corresponding new wavelengths
   
    :Returns:
    
    spectrum: array 
        interpolated to new wavelengths
    
    """        
    a,b,c = interpolate.splrep(wave, spec, k = 3)
    return interpolate.splev(new_wave,(a,b,c))


def calcStraightLine(x0, y0, x1, y1):
    """
    
    calculates the slope and the axis intercept of a straight line using to
    points in an x, y grid
    
    :Parameters:
    
    x0, y0, x1, y1: float
        2 points in an x, y grid
   
    :Returns:
    
    m, b: float 
        the slope and the axis intercept of the calculated straight line
    
    """
    m = (y1 - y0) /(x1 - x0) 
    b = (x1*y0 - x0*y1)/(x1 - x0)
    return m, b

def _get_weights( error, weighted, max_weight_gain ):
    """
    
    Helper function to convert an error array into a weights array.
    
    """
    if weighted and (error is not None):
        error = np.asarray(error)
        if error.min() > 0.0:
            # Prevent very tiny error estimates from skewing the result by
            # ignoring abnormally large weights.
            # TODO: Is this adjustment sensible? max_weight_gain is arbitrary.
            error_limit = error.mean() / max_weight_gain
            toosmall = np.where(error < error_limit)
            error[toosmall] = error_limit
            
            weights = 1.0 / error     
        else:
            weights = None      
    else:
        weights = None
    return weights

def _flatten_error_array( error, axis ):
    """
    
    Helper function to flatten an error array by taking the
    root mean square of the contributing pixels rather than
    the mean.
    
    """
    if error is not None:
        error_squared = error * error
        flat_errorsq = np.average( error_squared, axis=axis )
        flat_error = np.sqrt( flat_errorsq )
    else:
        flat_error = None
    return flat_error

def _flatten_quality_array( quality, axis, perfect ):
    """
    
    Helper function to flatten a quality array, with the result
    depending on whether perfection is required or not.
    
    """
    if quality is not None:
        if perfect:
            # Perfection is needed.
            # The overall quality is the maximum of the contributing quality.
            flat_quality = np.max( quality, axis=axis)
        else:
            # Perfection isn't needed.
            # The overall quality is the minimum of the contributing quality.
            flat_quality = np.min( quality, axis=axis)
    else:
        flat_quality = None
        
    return flat_quality

def spec2d_to_spec1d(data, error, quality, waveaxis=0, weighted=True,
                     max_weight_gain=10000.0, perfect=False):
    """
        
    Flatten a 2-D spectrum in the spatial direction to generate one
    combined spectrum.
        
    The output signal is the (weighted) average of the input signal.
    The output error is the RMS of the input error.
    The output quality is a combination of the input quality.
    
    :Parameters:
        
    data: array-like
        The 2-D spectrum to be flattened.
    error: array-like
        The 2-D error array associated with the data array.
    quality: array-like
        The 2-D data quality array associated with the data array.
        Good data elements are assumed to have quality=0.
    waveaxis: int, optional, default=0
        The wavelength axis for the data arrays.
    weighted: bool, optional
        Set to True to use the error array (if present) to make a
        weighted average.
        Set to False for an unweighted average.
    max_weight_gain: float, optional, default=10000.0
        Small errors can't weight a pixel by more than the average
        weighting multiplied by this amount.
    perfect: bool, optional
        Set to True to regard as bad any pixel in the spectrum
        which has contributing bad pixel.
        Set to False to regard a pixel in the spectrum as bad only
        if all of the contributing pixels are bad (and use the
        error array to mark the reduction in quality).

    :Returns:

    (flat_data, flat_error, flat_quality)
        A tuple containing a triplet of data arrays giving the signal,
        error and quality arrays for the extracted spectrum.
        
    """
    if waveaxis == 0:
        spatialaxis = 1
    else:
        spatialaxis = 0
        
    # If required, the signal average is weighted by the inverse of the error.
    weights = _get_weights( error, weighted, max_weight_gain )

    data = np.asarray(data)
    flat_data = np.average( data, axis=spatialaxis, weights=weights)
    flat_error = _flatten_error_array( error, axis=spatialaxis)
    flat_quality = _flatten_quality_array( quality, axis=spatialaxis,
                                           perfect=perfect)

    return (flat_data, flat_error, flat_quality)

def spec3d_to_spec1d(data, error, quality, waveaxis=0, weighted=True,
                     max_weight_gain=10000.0, perfect=False):
    """
        
    Flatten a spectral data cube in the spatial directions to
    generate one combined spectrum.
        
    The output signal is the (weighted) average of the input signal.
    The output error is the RMS of the input error.
        
    Not a particularly useful function unless the entire data cube
    contains the spectrum of one object.
    
    :Parameters:
        
    data: array-like
        The 2-D spectrum to be flattened.
    error: array-like
        The 2-D error array associated with the data array.
    quality: array-like
        The 2-D data quality array associated with the data array.
        Good data elements are assumed to have quality=0.
    waveaxis: int, optional, default=0
        The wavelength axis for the data arrays.
    weighted: bool, optional
        Set to True to use the error array (if present) to make a
        weighted average.
        Set to False for an unweighted average.
    max_weight_gain: float, optional, default=10000.0
        Small errors can't weight a pixel by more than the average
        weighting multiplied by this amount.
    perfect: bool, optional
        Set to True to regard as bad any pixel in the spectrum
        which has contributing bad pixel.
        Set to False to regard a pixel in the spectrum as bad only
        if all of the contributing pixels are bad (and use the
        error array to mark the reduction in quality).

    :Returns:

    (flat_data, flat_error, flat_quality)
        A tuple containing a triplet of data arrays giving the signal,
        error and quality arrays for the extracted spectrum.
        
    """
    # Determine a formula which averages every axis except the
    # wavelength axis.
    if waveaxis == 0:
        axis_list = (-1, -1)
    elif waveaxis == 1:
        axis_list = (-1, 0)
    else:
        axis_list = (1, 0)
        
    # If required, the signal average is weighted by the inverse of the error.
    weights1 = _get_weights( error, weighted, max_weight_gain )

    data = np.asarray(data)
    flat_data1 = np.average( data, axis=axis_list[0], weights=weights1)
    flat_error1 = _flatten_error_array( error, axis=axis_list[0])
    flat_quality1 = _flatten_quality_array( quality, axis=axis_list[0],
                                            perfect=perfect)

    weights2 = _get_weights( flat_error1, weighted, max_weight_gain )

    flat_data2 = np.average( flat_data1, axis=axis_list[1], weights=weights2)
    flat_error2 = _flatten_error_array( flat_error1, axis=axis_list[1])
    flat_quality2 = _flatten_quality_array( flat_quality1, axis=axis_list[1],
                                            perfect=perfect)
             
    return (flat_data2, flat_error2, flat_quality2)


def spec3d_to_image(data, error, quality, waveaxis=0, weighted=True,
                    max_weight_gain=10000.0, perfect=False):
    """
        
    Flatten a spectral data cube in the wavelength direction to
    generate one combined image.
        
    The output signal is the (weighted) average of the input signal.
    The output error is the RMS of the input error.
        
    Not a particularly useful function unless the entire data cube
    contains the spectrum of one object.
    
    :Parameters:
        
    data: array-like
        The 2-D spectrum to be flattened.
    error: array-like
        The 2-D error array associated with the data array.
    quality: array-like
        The 2-D data quality array associated with the data array.
        Good data elements are assumed to have quality=0.
    waveaxis: int, optional, default=0
        The wavelength axis for the data arrays.
    weighted: bool, optional
        Set to True to use the error array (if present) to make a
        weighted average.
        Set to False for an unweighted average.
    max_weight_gain: float, optional, default=10000.0
        Small errors can't weight a pixel by more than the average
        weighting multiplied by this amount.
    perfect: bool, optional
        Set to True to regard as bad any pixel in the spectrum
        which has contributing bad pixel.
        Set to False to regard a pixel in the spectrum as bad only
        if all of the contributing pixels are bad (and use the
        error array to mark the reduction in quality).

    :Returns:

    (flat_data, flat_error, flat_quality)
        A tuple containing a triplet of data arrays giving the signal,
        error and quality arrays for the flattened image.
        
    """
    # If required, the signal average is weighted by the inverse of the error.
    weights = _get_weights( error, weighted, max_weight_gain )

    data = np.asarray(data)
    flat_data = np.average( data, axis=waveaxis, weights=weights)
    flat_error = _flatten_error_array( error, axis=waveaxis) 
    flat_quality = _flatten_quality_array( quality, axis=waveaxis,
                                           perfect=perfect)

    return (flat_data, flat_error, flat_quality)



# A minimal test is run when this file is run as a main program.
if __name__ == '__main__':
    print("Testing the module\n")
    PLOTTING = True     # Set to False to turn off plotting.

    print("condense")
    spec = np.arange(50, 100, 0.2)
    new_spec = condense(spec, 5)
    assert len(spec) == len(new_spec)
    print(new_spec)
    new_spec = condense(spec, 10, method = "Mean",running = False)
    print(new_spec)
    assert len(spec)/10 == len(new_spec)
    
    print("convGauss")
    x = np.linspace(0, 49)
    y = convGauss(x, 5, 1)
    assert np.mean(x) == np.mean(y)
    assert len(y) == len(x) - 4

    # Construct a simulated 7 x 9 array in which wavelength runs in the
    # Y direction (waveaxis=0 or -2)
    waveaxis = 0
    a = [[10, 10, 20, 40,  20, 10, 10],
         [10, 10, 20, 50,  20, 10, 10],
         [10, 10, 20, 150,  20, 10, 10],
         [10, 10, 20, 70,  50, 25, 10],
         [10, 10, 20, 50,  20, 10, 10],
         [10, 10, 40, 70,  30, 10, 10],
         [10, 10, 20, 170,  20, 10, 10],
         [10, 10, 20, 60,  20, 10, 10],
         [10, 10, 20, 40,  20, 10, 10]]
    b = [[0.5, 0.5, 0.2, 0.05, 0.2, 0.1, 0.5],
         [0.6, 0.5, 0.2, 0.05, 0.2, 0.5, 0.6],
         [0.6, 0.5, 0.2, 0.05, 0.2, 0.5, 0.6],
         [0.5, 0.4, 0.2, 0.05, 0.2, 0.4, 0.5],
         [0.4, 0.5, 0.2, 0.02, 0.2, 0.5, 0.5],
         [0.1, 0.4, 0.2, 0.02, 0.2, 0.6, 0.4],
         [0.6, 0.5, 0.2, 0.03, 0.2, 0.5, 0.6],
         [0.6, 0.5, 0.2, 0.02, 0.2, 0.5, 0.5],
         [0.5, 0.5, 0.2, 0.05, 0.2, 0.5, 0.4]]
    
    print("Creating slope product from arrays")
    product = MiriMeasuredModel( data=a, err=b )
    print(product)
    if PLOTTING:
        product.plot()
    
    print("lrs2D_spextract")
    lrs2d = lrs2D_spextract(product, xmin = 3, xmax = 6, ymin = 3, ymax = 6)
    print(lrs2d)
    if PLOTTING:
        lrs2d.plot()
    
    print("subtractBackground")   
    lrs2d = subtractBackground(lrs2d, lrs2d)
    print(lrs2d)
    if PLOTTING:
        lrs2d.plot()
    
    print("lrs_extract_spec")
    spec = lrs_extract_spec(lrs2d)
    print(spec)
    if PLOTTING:  
        spec.plot()    
    print("optimalSpecExtraction")
    sp = optimalSpecExtraction(product.copy())
    print(sp)
    if PLOTTING:  
        sp.plot()    
    print("lrs_extract_spec_with_fit")
    sp = lrs_extract_spec_with_fit(product.copy(), verbose = True, makeplot=True)
    print(sp)         
    
    print("get_psf_fit")
    sp = get_psf_fit(product.copy(), verbose = True, makeplot=True)
    print(sp)         
    
    
    del product, lrs2d, spec , sp
    
    print("interpolWaveOnRows")
    wave = np.arange(5, 10, 0.5)
    print("original wavelengths ", wave)
    pos = np.arange(0, 5, 0.5)
    rows, new_wave = interpolWaveOnRows(pos, wave)
    print("rows ", rows)
    print("new wavelengths ",new_wave)
    
    print("interpolLin")
    spec = np.arange(100, 105, 0.5)
    print("original spectrum ", spec)
    new_spec = interpolLin(wave, spec, new_wave)
    print("new spectrum ", new_spec)
    
    print("interpol_lin")
    spec = np.arange(100, 105, 0.5)
    print("original spectrum ", spec)
    new_spec = interpol_lin(wave, spec, new_wave)
    print("new spectrum ", new_spec)            
    
    print("interpolSpline")
    spec = np.arange(100, 105, 0.5)
    print("original spectrum ", spec)
    new_spec = interpolSpline(wave, spec, new_wave)
    print("new spectrum ", new_spec)                
   
    print("calcStraightLine")
    print(calcStraightLine(1,2,3,4))
    print("Test finished.")

