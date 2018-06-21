#!/usr/bin/env python
#
# Script 'convert_im_flux_to_phot'
#
# :History:
# 
# 05 Nov 2015: Created
#
# @author: Steven Beard (UKATC)
#
"""

Script `convert_im_flux_to_phot` converts a data model in the old
MiriImagingFluxconversionModel format into a data model in the new
MiriImagingPhotometricModel format. The main differences are::

    * The old data model contains a table with 3 columns::

        FILTER - The name of an imager filter
        FACTOR - The flux conversion factor for this filter, in
            (Jy/arcsec2)/(DN/s/pixel)) units.
        UNCERTAINTY - The uncertainty in the flux conversion factor.
        
    * The new data model contains a table with 7 columns::
    
        FILTER - The name of an imager filter
        SUBARRAY - The name of a subarray mode (empty string if the flux
            conversion factor does not depend on subarray mode).
        PHOTMJSR - The flux conversion factor for this filter/subarray
            combination, in (MJ/sr)/(DN/s/pixel) units.
        UNCERTAINTY - The uncertainty in the flux conversion factor.
        NELM - Number of elements in WAVELENGTH array.
        WAVELENGTH - Wavelength array.
        RELRESPONSE - Array containing relative flux conversion factors
            as a function of wavelength
        
        For MIRI imager photometric models, the NELM, WAVELENGTH and
        RELESPONSE columns are all set to zero.

    * The new data model also contains additional items of metadata::
    
        PIXAR_SR - The nominal pixel area for the detector in steradians.
        PIXAR_A2 - The nominal pixel area for the detector in square arcseconds


The following command arguments are defined by position:

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out" appended.

The command also takes the following options:

    --subarray <subarray-name-string>:
        The subarray to which the conversion factors apply.
        Defaults to no subarray.
    --factor <conversion-factor>:
        The multiplication factor needed to convert (Jy/arcsec2) into (MJ/sr).
        Defaults to the constant given below.
    --pixar_a2 <pixel-area>:
        The pixel area to write to the PIXAR_A2 keyword.
        Defaults to the constant given below.
    --verbose or -v:
        Generate more output.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""



import optparse
import os, sys, time

#import astropy.io.ascii as ascii
import numpy as np

from miri.datamodels.miri_fluxconversion_models import \
    MiriImagingFluxconversionModel
from miri.datamodels.miri_photometric_models import \
    MiriImagingPhotometricModel, ARCSEC2_PER_STERADIAN

def convert_model( input_model, subarray, factor, pixar_a2, pixar_sr):
    """
    
    Converts the data model
    
    """
    assert isinstance(input_model, MiriImagingFluxconversionModel)
    
    # Extract the FILTER, FACTOR and UNCERTAINTY columns from the old
    # flux_table, create a new flux table with the flux column multiplied
    # by the given conversion factor.
    flux_table = input_model.flux_table
    new_flux_table = []
    for element in flux_table:
        filter = element[0]
        subarray = subarray
        photmjsr = element[1] * factor
        uncertainty = element[2] * factor
        new_flux_table.append( (filter, subarray, photmjsr, uncertainty) )
        
    # Create a data model in the new format:
    output_model = MiriImagingPhotometricModel(phot_table=new_flux_table,
                                               pixar_a2=pixar_a2,
                                               pixar_sr=pixar_sr)
    # Copy the metadata from the original data model, but ensure the
    # output data type is preserved.
    output_model.copy_metadata(input_model)
    output_model.meta.reftype = 'PHOTOM'
    return output_model

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts an old MIRI imager flux conversion model into a "
    usage += "standard MIRI CDP-5 format MIRI imager photometric model."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--subarray", dest="subarray", type="str", action="store",
                      default="", help="Subarray to which the flux data applies."
                     )
    deffactor = ARCSEC2_PER_STERADIAN / 1.0e6
    help_text += \
        "The default flux conversion factor (--factor) is %f\n" % \
        deffactor
    parser.add_option("", "--factor", dest="factor", type="float", action="store",
                      default=deffactor, help="Multiplication factor for flux units."
                     )
    defarea = 0.11 * 0.11
    help_text += \
        "The default pixel area (--pixar_a2) is %f square arcseconds.\n" % \
        defarea
    parser.add_option("", "--pixar_a2", dest="pixar_a2", type="float", action="store",
                      default=defarea, help="Pixel area in arcseconds squared."
                     )
#     parser.add_option("-p", "--plot", dest="makeplot", action="store_true",
#                       help="Plot flux table"
#                      )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the FITS file if it already exists"
                     )
    (options, args) = parser.parse_args()

    try:
        inputfile = args[0]
        if len(args) > 1:
            outputfile = args[1]
        else:
            outputfile = inputfile + "_out.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    subarray = options.subarray
    factor = options.factor
    pixar_a2 = options.pixar_a2
    pixar_sr = pixar_a2 / ARCSEC2_PER_STERADIAN
    verb = options.verb
#     makeplot = options.makeplot
    overwrite = options.overwrite

    # Read the input data model.
    if verb:
        print("Reading %s" % inputfile)

    with MiriImagingFluxconversionModel( inputfile ) as input_model:
        if verb:
            print("Input data model")
            print(input_model)
#         if makeplot:
#             input_model.plot()
        if verb:
            strg = "Converting model "
            if subarray:
                strg += "for subarray %s " % subarray
            strg += "with flux column multiplier=%f, " % factor
            strg += "PIXAR_A2=%f and PIXAR_SR=%f." % (pixar_a2, pixar_sr)
            print(strg)
        output_model = convert_model( input_model, subarray, factor,
                                      pixar_a2, pixar_sr)
        # Update the reference file HISTORY and version control records
        # within the new model.
        output_model.meta.description = \
            'CDP-5 MIRI Imager Flux Calibration Factors'
        output_model.meta.version = '05.00.00'
        output_model.meta.useafter = '2015-11-20'
        output_model.add_history("DOCUMENT: MIRI-TR-00009-UoC Issue2.2.pdf")
        if verb:
            print("Output data model")
            print(output_model)

        output_model.save( outputfile, overwrite=overwrite)
        print("Data saved to %s\n" % outputfile)
        del input_model
        del output_model
