#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# 27 May 2015: Replaced pyfits with astropy.io.fits


"""

Checking and debugging the STScI data model (jwst_lib.models).

This module demonstrates problems found while using the STScI data model
and shows the solutions to those problems recommended by the STScI team.

The module can be rerun to make sure the solutions still work.

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import numpy as np
#import numpy.ma as ma
print("Using numpy V" + str(np.__version__))

import astropy.io.fits as pyfits
print("Using astropy.io.fits (pyfits) V" + str(pyfits.__version__))

# Import the STScI image model (and utilities if required)
#import jwst_lib.models.util as jmutil
#import jwst_lib.models.schema as mschema
#from jwst_lib.models.model_base import DataModel
from jwst_lib.models.image import ImageModel
from jwst_lib.models.flux import FluxModel
#from miri.datamodels.miri_fluxconversion_models import FluxModel
# from jwst_lib.models.ramp import RampModel

# from miri.datamodels.miri_fluxconversion_model import \
#     MiriFluxconversionModel

from miri.datamodels.miri_badpixel_model import MiriBadPixelMaskModel
from miri.datamodels.miri_measured_model import MiriMeasuredModel
from miri.datamodels.miri_fluxconversion_models import MiriImagingFluxconversionModel

def save_old_data( filename, data, header=None, overwrite=False, add_comment=False ):
    """
    
    Save data to an old format FITS file with data contained in
    the PRIMARY HDU.
    
    """
    fitsheader = pyfits.Header()
    fitsheader.update("FILENAME", filename, "Name of file")
    
    # Add a COMMENT card to the primary header, if requested
    if add_comment:
        fitsheader.add_comment("This is a comment in the original header.")
            
    if header is not None:
        for key in header.keys():
            fitsheader.update(key, header[key], "")
            
    hdulist = pyfits.HDUList([])
    primary_hdu = pyfits.PrimaryHDU(data=data, header=fitsheader)
    hdulist.append(primary_hdu)
    hdulist.writeto(filename, overwrite=overwrite)
    del hdulist
            
def load_old_data( filename ):
    """
    
    Reads a header and data from an old format file where
    data are contained in the PRIMARY HDU.
    
    """
    # Read the FITS header and primary data array from the FITS file.
    try:
        hdulist = pyfits.open(filename)
        fitsheader = hdulist[0].header
        data = hdulist[0].data
        
    finally:
        hdulist.close()
        del hdulist

    return fitsheader, data

#
# The following code is executed only when this file is run as a main program.
#
if __name__ == '__main__':
    print("Checking the STScI data model module (jwst_lib.models).")
    
    VERBOSE = False
    
    FILE_LOCKING_PROBLEM = False
    FILE_LOCKING_SOLUTION1 = True
    DATA_CONVERSION_PROBLEM = False
    DATA_CONVERSION_SOLUTION = True
    WORLD_COORDINATES_PROBLEM = False
    TABLE_ARRAY_PROBLEM = True

    # -----------------------------------------------------------------------
    #
    # Demonstrates the file locking problem.
    #
    if FILE_LOCKING_PROBLEM:
        data = np.array([[1,2,3],[4,5,6],[7,8,9]])
        filename = "testfile.fits"
        model = ImageModel( data=data )
        model.save(filename, overwrite=True)
        readback = ImageModel(filename)
        assert model.data.all() == readback.data.all()
        del readback, model
        # This line fails because the file is still locked.
        os.remove(filename)
    #
    # Problem initially solved by using the "with" context manager.
    # (However, the problem returned in April 2013 and might need
    # a different solution.)
    #
    if FILE_LOCKING_SOLUTION1:
        data = np.array([[1,2,3],[4,5,6],[7,8,9]])
        filename = "testfile.fits"
        with ImageModel( data=data ) as model:
            model.save(filename, overwrite=True)
            with ImageModel(filename) as readback:
                assert model.data.all() == readback.data.all()
                del readback
            del model
        # This line used to work.
        os.remove(filename)

    # -----------------------------------------------------------------------
    #
    # Demonstrates the data model conversion problem.
    #
    infile = "test_input_file.fits"
    outfile = "test_output_file.fits"
        
    # Create a temporary file in old format, with data contained
    # in the PRIMARY HDU
    data2D = [[1,2,3],[4,5,6],[7,8,9]]
    image_data = np.array(data2D)
    save_old_data( infile, image_data, overwrite=True)
            
    # Read back the temporary file created above
    header, data = load_old_data( infile )
    
    if DATA_CONVERSION_PROBLEM:
        # Attempt to convert that file into a new mask object.
        # The original file is read as a template to preserve all the
        # existing metadata in the PRIMARY HDU.
        with MiriBadPixelMaskModel( init=infile ) as newproduct:
            # The mask array is explicitly defined using the
            # array read from the PRIMARY HDU.
            newproduct.mask = data
            # The print statement shows the data product looks ok.
            if VERBOSE:
                print("\nContents of the data product created from", infile)
                print((30 * "-x") + "-")
                print(str(newproduct))
            # But an attempt to save it to a file results in the exception
            # "AttributeError: 'numpy.ndarray' object has no attribute '_scale_back'"
            newproduct.save( outfile, overwrite=True )
            del newproduct

    if DATA_CONVERSION_SOLUTION:
        # The solution is to make a copy of the data product before
        # attempting to alter its contents.
        with MiriBadPixelMaskModel( init=infile ) as oldproduct:
            newproduct = oldproduct.copy()
            # The mask array is explicitly defined using the
            # array read from the PRIMARY HDU.
            newproduct.mask = data
            # The print statement shows the data product looks ok.
            if VERBOSE:
                print("\nContents of the data product cloned from", infile)
                print((30 * "-x") + "-")
                print(str(newproduct))
            # and this time the data product saves successfully.
            newproduct.save( outfile, overwrite=True )
            del newproduct
        os.remove(outfile)
    # Delete these objects to avoid the file locking problem with infile.
    del header,data
    os.remove(infile)

    # -----------------------------------------------------------------------
    #
    # Demonstrates the world coordinates conversion problem.
    #
    infile = "test_input_file.fits"
    outfile = "test_output_file.fits"
        
    # Create a temporary file in old format, with data contained
    # in the PRIMARY HDU. This time give the primary header some
    # world coordinates keywords.
    data2D = [[1,2,3],[4,5,6],[7,8,9]]
    image_data = np.array(data2D)
     
    header = {'CRPIX1' : 1,
              'CRVAL1' : 0.0,
              'CDELT1' : 0.5,
              'CTYPE1' : 'alpha',
              'CUNIT1' : 'arcsec',
              'CRPIX2' : 1,
              'CRVAL2' : 0.0,
              'CDELT2' : 0.5,
              'CTYPE2' : 'beta',
              'CUNIT2' : 'arcsec'
              }
     
    save_old_data( infile, image_data, header=header, overwrite=True, add_comment=True)
             
    # Read back the temporary file created above
    header, data = load_old_data( infile )
     
    if WORLD_COORDINATES_PROBLEM:
        # Attempt to convert that file into a new measured data object.
        # The original file is read as a template to preserve all the
        # existing metadata in the PRIMARY HDU.
        with MiriMeasuredModel( init=infile ) as oldproduct:
            newproduct = oldproduct.copy()
            # The data array is explicitly defined using the
            # array read from the PRIMARY HDU.
            newproduct.data = data
            # The print statement shows the data product looks ok, except
            # the world coordinates keywords have been moved into the
            # "_extra_fits" location for the PRIMARY HDY.
            if VERBOSE:
                print("\nContents of the data product with WCS cloned from", infile)
                print((30 * "-x") + "-")
                print(str(newproduct))
            # Saving the data product to FITS at this point would leave
            # the world coordinates in the wrong HDU.

            # Copy the WCS info from the PRIMARY HDU to the SCI HDU.
            # This line can fail with
            # "TypeError: __init__() got an unexpected keyword argument 'fix'"
            # until the fix=True argument is removed from jwst_lib.wcs.
            # The solution to the problem is to use the correct version
            # of pywcs.
            # It can also fail with:
            # RuntimeError: module compiled against API version 7 but this version of numpy is 6
            wcs = newproduct.get_fits_wcs( hdu_name='PRIMARY' )
            newproduct.set_fits_wcs( wcs, hdu_name='SCI')

            # The print statement now shows the WCS keywords have been
            # written to the _extra_fits area for the SCI HDU, as expected.
            if VERBOSE:
                print("\nContents of the data product after copying the WCS")
                print((30 * "-x") + "-")
                print(str(newproduct))
            # The data product can now be saved with the correct WCS info.
            newproduct.save( outfile, overwrite=True )
            del newproduct
            
        os.remove(outfile)
#       Delete these objects to avoid the file locking problem with infile.
        del header,data
        os.remove(infile)
        
    # -----------------------------------------------------------------------
    #
    # Demonstrates the table data problems.
    #
    testfile = "test_table_data.fits"
    testfile2 = "test_table_data_copy.fits"

    flux_im = [('F560W',  1.0e-5,  1.0e-7),
        ('F770W',  1.1e-5,  1.6e-7),
        ('F1000W', 1.2e-5,  1.1e-7),
        ('F1130W', 1.3e-5,  0.0),
        ('F1280W', 1.4e-5,  1.7e-7),
        ('F1500W', 1.5e-5,  2.0e-7),
        ('F1800W', 1.6e-5,  1.1e-7),
        ('F2100W', 1.7e-5,  3.0e-7),
        ('F2550W', 1.8e-5,  4.2e-7),
        ('F2550WR',1.8e-5,  4.2e-7),
        ('F1065C', 2.8e-5,  4.2e-7),
        ('F1140C', 6.8e-6,  8.2e-8),
        ('F1550C', 2.1e-5,  2.2e-7),
        ('F2300C', 9.8e-6,  5e-8),
        ]
  
    if TABLE_ARRAY_PROBLEM:            
        # Create a standard data product containing table data.
        print("Creating flux model and saving to file", testfile)
        with FluxModel( flux_table=flux_im ) as datamodel:
            print(datamodel.flux_table)
            datamodel.save(testfile, overwrite=True)    
            del datamodel
            
        # Read the model back from the file.
        print("Reading flux model back from the same file and copying to",
              testfile2)
        with FluxModel( init=testfile ) as newmodel:
            # The following print statement triggers the problem.
            # The solution is to use a recent version of pyfits
            if str(pyfits.__version__) >= '3.1.2':
                print(newmodel.flux_table)
            # This line fails with "AttributeError: 'numpy.ndarray' object
            # has no attribute '_scale_back'"
            newmodel.save(testfile2, overwrite=True)
            del newmodel        

    # -----------------------------------------------------------------------
    #
    # Demonstrates the table problem with old data.
    # Modify the input file name as required.
    #
    testfile = "D:/MIRI_FM_IM_AbsFlux_01.00.00.fits"
    testfile2 = "flux_table_data_copy.fits"
    testfile3 = "flux_table_data_copy2.fits"
    if TABLE_ARRAY_PROBLEM and os.path.isfile(testfile):            
        # Read a flux model back from an old existing file.
        print("Reading flux model back from file", testfile)
        with MiriImagingFluxconversionModel( init=testfile ) as fluxmodel:
            # The following print statement triggers the problem.
            # The solution is to use a recent version of pyfits
            if str(pyfits.__version__) >= '3.1.2':
                print(fluxmodel.flux_table)
            # This line fails with "AttributeError: 'numpy.ndarray' object
            # has no attribute '_scale_back'"
            fluxmodel.save(testfile2, overwrite=True)
            del fluxmodel        

        # Read the model back from the file just created.
        print("Reading flux model back from the same file and copying to",
              testfile3)
        with MiriImagingFluxconversionModel( init=testfile2 ) as newmodel:
            # The following print statement triggers the problem.
            # The solution is to use a recent version of pyfits
            if str(pyfits.__version__) >= '3.1.2':
                print(newmodel.flux_table)
            # This line fails with "AttributeError: 'numpy.ndarray' object
            # has no attribute '_scale_back'"
            newmodel.save(testfile3, overwrite=True)
            del newmodel        

    print("Checks complete.")
    
