#!/usr/bin/env python
#
# :History:
# 
# 04 Oct 2018: Created
# 09 Oct 2018: Convert the PHASE2 and PHASE3 data.
# 10 Oct 2018: Added --version and --document as parameters. Delete the
#              resolving_power attribute left over from the old data model.
# 17 Oct 2018: 'ANY' replaced by 'N/A' in metadata wildcards.
# 24 Oct 2018: Add the missing TELESCOP and INSTRUME keywords
#
# @author: Steven Beard (UKATC)
#

"""

Script `convert_mrs_resolution` creates a FITS file or MRS spectral resolution
data in standard MIRI/JWST data model format from a FITS file in Jeb Bailey's
original format.

The script relies on the MIRI MRS spectral resolution data product, which is
based on the STScI JWST data models.

The following command arguments are defined by position:

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out" appended.

The command also takes the following options:

    --document
        Override the default document name and provide a new one.
    --version
        Override the version number in the file and provide a new one.
    --verbose or -v:
        Generate more output.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""

import optparse
import os, sys, time
import warnings

import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.miri_spectral_spatial_resolution_model \
    import MiriMrsResolutionModel, MAX_NELEM

def wildcard_filter( input_string ):
    """
    
    A helper function which filters out a wildcard string containing
    'ANY' and converts it to 'N/A'.
    
    """
    if str(input_string).strip() == 'ANY':
        return 'N/A'
    else:
        return input_string

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile [outputfile]\n"
    usage += "Converts a Jeb Bailey format MRS spectral resolution file into "
    usage += "standard MIRI CDP-7 format."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--version", dest="version", type="string",
                      default=None, help="New version number (overriding file contents)"
                     )
    parser.add_option("", "--document", dest="document", type="string",
                      default=None, help="New document name (overriding default)"
                     )
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

    document = options.document
    version = options.version
    verb = options.verb
    overwrite = options.overwrite

    # Read the header and data from given file.
    if verb:
        strg = "Reading %s" % inputfile
        print(strg)
    
    with pyfits.open(inputfile) as hdulist:

        # Read the primary header
        fitsheader = hdulist[0].header
        
        # Read the spectral resolution tables and coefficients.
        try:
            psf_fwhm_alpha = hdulist['PSF_FWHM_ALPHA'].data
            psf_fwhm_beta = hdulist['PSF_FWHM_BETA'].data
        except KeyError as e:
            strg = "%s does not contain MIRI MRS spectral resolution data." % inputfile
            strg += "\n  %s" % str(e)
            raise TypeError(strg)

        try:
            resol_data = hdulist['RESOL_DATA'].data
            mlsf_data = hdulist['MLSF_DATA'].data
            phase1_data = hdulist['PHASE1_DATA'].data
            phase2_data_orig = hdulist['PHASE2_DATA'].data
            phase3_data_orig = hdulist['PHASE3_DATA'].data
            etalon_data = hdulist['ETALON_DATA'].data
        except KeyError as e:
            strg = "%s contains pre-CDP-7 MIRI MRS spectral resolution data." % inputfile
            strg += "\n  Please give the name of a Jeb Bailey MRS spectral resolution file."
            strg += "\n  %s" % str(e)
            raise TypeError(strg)
        
        # The PHASE2 and PHASE3 data tables must be converted to the new format.
        phase2_data = []
        for phase2_row in phase2_data_orig:
            if len(phase2_row) <= 4:
                strg = "%s is not a Jeb Bailey MRS spectral resolution file." % inputfile
                strg += "\n  PHASE2_DATA has only %d columns. " % len(phase2_row)
                strg += "Has this file already been converted?"
                raise TypeError(strg)
            phase2_domain_low = phase2_row[0]
            phase2_domain_high = phase2_row[1]
            phase2_norder = len(phase2_row) - 2
            phase2_newcoeff = MAX_NELEM * [0.0]
            for coeff in range(0, phase2_norder):
                phase2_newcoeff[coeff] = phase2_row[2+coeff]
            phase2_data.append( (phase2_domain_low, phase2_domain_high, phase2_norder, phase2_newcoeff) )
        
        phase3_data = []
        for phase3_row in phase3_data_orig:
            phase3_norder = len(phase3_row)
            phase3_newcoeff = MAX_NELEM * [0.0]
            for coeff in range(0, phase3_norder):
                phase3_newcoeff[coeff] = phase3_row[coeff]
            phase3_data.append( (phase3_norder, phase3_newcoeff) )

        # Copy the contents of this file to a new data model
        with MiriMrsResolutionModel( psf_fwhm_alpha=psf_fwhm_alpha,
                                     psf_fwhm_beta=psf_fwhm_beta,
                                     resol_data=resol_data, mlsf_data=mlsf_data,
                                     phase1_data=phase1_data,
                                     phase2_data=phase2_data,
                                     phase3_data=phase3_data,
                                     etalon_data=etalon_data ) as resolmodel:

            # Copy over the metadata.
            origin = 'MIRI European Consortium'
            author = fitsheader['AUTHOR']
            pedigree = fitsheader['PEDIGREE']
            if version is None or not version:
                version = fitsheader['VERSION']
            date = fitsheader['DATE']
            useafter = fitsheader['USEAFTER']
            description = fitsheader['DESCRIP']
            resolmodel.set_housekeeping_metadata(origin, author=author,
                                pedigree=pedigree, version=version, date=date,
                                useafter=useafter, description=description )
            resolmodel.meta.filename_original = os.path.basename(inputfile)
            resolmodel.meta.filename = os.path.basename(outputfile)
            
            modelnam = wildcard_filter( fitsheader['MODELNAM'] )
            detector = wildcard_filter( fitsheader['DETECTOR'] )
            detsetng = wildcard_filter( fitsheader['DETSETNG'] )
            if 'BAND' in fitsheader:
                band = wildcard_filter( fitsheader['BAND'] )
            else:
                band = 'N/A'
            if 'CHANNEL' in fitsheader:
                channel = wildcard_filter( fitsheader['CHANNEL'] )
            else:
                channel = 'N/A'
            filt = 'N/A'
            resolmodel.set_telescope()
            resolmodel.set_instrument_metadata(detector, modelnam=modelnam,
                        detsetng=detsetng, filt=filt, channel=channel,
                        band=band)

            resolmodel.meta.exposure.readpatt = wildcard_filter( fitsheader['READPATT'] )
            resolmodel.meta.exposure.nframes = 1
    
            if 'SUBARRAY' in fitsheader:
                subarray = fitsheader['SUBARRAY']
            else:
                subarray = 'GENERIC'
            resolmodel.set_subarray_metadata( subarray )
            resolmodel.set_exposure_type( detector=detector, subarray=subarray )
            
            if document is None or not document:
                doc_strg = 'DOCUMENT: MIRI-RP-00514-NLC FM-MRS-Spectral-Resolution Issue 9/3/18'
            else:
                doc_strg = 'DOCUMENT: %s' % document
            resolmodel.add_history(doc_strg)
            resolmodel.add_history('SOFTWARE: IDL and Python')
            resolmodel.add_history('DATA USED: Derived from FM data')
            resolmodel.add_history('DIFFERENCES: New data model format')
            
            # Delete the RESOLVING_POWER HDU left over from the old data model
            # before saving the file.
            if hasattr(resolmodel, "resolving_power"):
                del resolmodel.resolving_power
            if verb:
                print(resolmodel)
    
            resolmodel.save( outputfile, overwrite=overwrite)
            if verb:
                print("Data saved to %s\n" % outputfile)
        del resolmodel
    hdulist.close()
    del hdulist
