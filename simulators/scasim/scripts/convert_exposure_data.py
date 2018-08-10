#!/usr/bin/env python
#
# Script 'convert_exposure_data' converts between the two kinds of
# exposure data supported by SCASim. There are two variants:
#
# level1b_to_dhas: Converts STScI "level 1b" ramp data into the legacy
#                  FITSWriter format used by DHAS.
#
# dhas_to_level1b: Converts the FITSWriter format used by DHAS into
#                  STScI "level 1b" ramp data.
#
# :Caveats:
#
# FITSWriter format does not support data quality or uncertainty arrays, so
# this information will not be propagated.
#
# The two data formats use different metadata standards, so the metadata
# conversion may not be entirely correct.
#
# :History:
#
# 23 Feb 2016: Original version, based on existing SCASim conversion utilities.
# 09 Mar 2016: Corrected typo in doc string.
# 14 Mar 2016: Convert DATE metadata fields explicitly to a string when
#              converting from ramp data to DHAS format.
# 04 May 2016: Modified to process the reference output data.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
# 10 Feb 2017: Corrected typo --> level2 should be level1.
# 07 Aug 2017: Explicitly translate some important DHAS metadata.
#              Corrected a bug: DETROWS keyword is not normally found in
#              level 1 ramp data. Ensure all attempts to read a FITS
#              keyword are within a try/except block.
# 04 Dec 2017: Added missing "overwrite" parameter to functions.
# 27 Apr 2018: Corrected exception syntax for Python 3.
# 17 May 2018: Python 3: Converted dictionary keys return into a list.
# 18 May 2018: Changed deprecated logger.warn() to logger.warning().
# 
# @author: Steven Beard (UKATC)
#

"""

Script 'convert_exposure_data' converts between the two kinds of
exposure data supported by SCASim.

NOTE: The metadata may not convert correctly. Check for warnings.

The compulsory command arguments are:

    inputfile
        The path+name of the input file.
    outputfile
        The path+name of the output file.

The following optional parameters may be provided by keyword:

    --convert
        The kind of conversion required.
        
        'level1b_to_dhas' converts STScI "level 1b" ramp data into
        the legacy FITSWriter format used by DHAS.
        
        'dhas_to_level1b' converts the FITSWriter format used by
        DHAS into STScI "level 1b" ramp data.
   
The command also takes the following options:

    --silent or -s:
        Generate no output.
    --verbose or -v:
        Generate more output and some data plots.
    --debug or -d:
        Generate maximum output.
    --overwrite or -o:
        Overwrite any existing SCA format FITS file.

"""

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.simulators.convert_exposure_data") # Get a default parent logger

import optparse
import os, sys, time
import re
import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels import MiriRampModel
from miri.datamodels.sim import MiriExposureModel
from miri.simulators.scasim.exposure_data import ExposureData, load_exposure_data

def convert_dhas_to_level1b( inputfile, outputfile, verbose=0, overwrite=False ):
    """
    
    Convert a DHAS/FitsWriter format file into a STScI level 1b
    ramp data file.
    
    """
    strg = "NOTE: Full metadata conversion is not guaranteed - check for warnings"
    LOGGER.info(strg)
    input_model = load_exposure_data( inputfile )
    if verbose > 1:
        print( input_model )

    # Determine the size of the reference output array.
    rows = input_model.data.shape[-2]
    columns = input_model.data.shape[-1]
    try:
        detrows = input_model.get_fits_keyword('ROWSTOP')
        toprows = rows - detrows
        if toprows > 0:
            refcolumns = toprows * columns / detrows
            refrows = detrows
        else:
        # No reference output
            refcolumns = 0
            refrows = 0
    except KeyError:
        # No reference output
        LOGGER.warning("No reference output data. Cannot deduce size of reference image.")
        refcolumns = 0
        refrows = 0
    
    output_model = MiriExposureModel( rows=input_model.rows,
                                      columns=input_model.columns,
                                      refrows=refrows, refcolumns=refcolumns,
                                      ngroups=input_model.ngroups,
                                      nints=input_model.nints,
                                      readpatt=input_model.readpatt )
    output_model.set_exposure(input_model.data, dq=None)
    
    # Copy the FITS header from the input model to the metadata
    # of the output model.
    # NOTE: No attempt is made to translate between different metadata standards.
    ignore = ['SIMPLE', 'EXTEND', 'BITPIX', 'NAXIS', 'BZERO', 'BSCALE', 'HISTORY', 'COMMENT']
    for fitskw in list(input_model._primary_header.keys()):
        addkw = True
        if not fitskw:
            continue
        for teststrg in ignore:
            if teststrg in fitskw:
                addkw = False
        if addkw:
            value = input_model._primary_header[fitskw]
            if fitskw == 'READOUT':
                # Translate the READOUT keyword into READPATT
                output_model.set_fits_keyword('READPATT', value)
            elif fitskw == 'MIRMODEL':
                # Translate the MIRMODEL keyword into MODELNAM
                if value.startswith('F') or value.startswith('f'):
                    value = 'FM'
                elif value.startswith('V') or value.startswith('v'):
                    value = 'VM'
                elif value.startswith('J') or value.startswith('j'):
                    value = 'JPL'
                else:
                    value = 'FM'
                output_model.set_fits_keyword('MODELNAM', value)
            elif fitskw == 'SCA_ID':
                # Special case. SCA_ID must be converted to string.
                output_model.set_fits_keyword('SCA_ID', str(value))
            elif fitskw == 'FRMRSETS':
                # Special case. FRMRSETS must be converted to integer.
                output_model.set_fits_keyword('FRMRSETS', int(value))
            elif fitskw == 'SUBARRAY':
                # Special case. SUBARRAY must be converted to string
                # and 'F' or 'false' are not allowed.
                value = str(value)
                if value.startswith('F') or value.startswith('f'):
                    value = 'FULL'
                output_model.set_fits_keyword('SUBARRAY', value)
            else:
                # Otherwise copy the keyword
                try:
                    output_model.set_fits_keyword(fitskw, value)
                except KeyError:
                    strg = "\n\tFITS keyword \'%s\' not found in data model. " % fitskw
                    LOGGER.warning(strg)
                except ValueError:
                    strg = "\n\tFailed to convert metadata item with "
                    strg += "FITS keyword \'%s\'. " % fitskw
                    LOGGER.error(strg)
                except Exception as e:
                    strg = "\n\tFailed to convert metadata item with "
                    strg += "FITS keyword \'%s\'. " % fitskw
                    strg += "\n\t%s" % str(e)
                    LOGGER.error(strg)
    
    hstrg = "Created from DHAS file %s" % inputfile
    output_model.add_history(hstrg)
    hstrg = "Converted by SCASim convert_exposure_data script."
    output_model.add_history(hstrg)
    output_model.save(outputfile, overwrite=overwrite)

def convert_level1b_to_dhas( inputfile, outputfile, verbose=0, overwrite=False ):
    """
    
    Convert a STScI level 1b ramp data file into a DHAS/FitsWriter
    format file.
    
    """
    strg = "NOTE: Full metadata conversion is not guaranteed - check for warnings"
    LOGGER.info(strg)
    input_model = MiriRampModel( inputfile )
    if verbose > 1:
        print( input_model )
    
    # Deduce the number of columns, rows, groups and integrations from
    # the primary data shape.   
    inputdata = input_model.data
    if inputdata.ndim >= 4:
        nints = inputdata.shape[0]
        ngroups = inputdata.shape[1]
        rows = inputdata.shape[2]
        columns = inputdata.shape[3]
    elif input_model.data.ndim >= 3:
        nints = 1
        ngroups = inputdata.shape[0]
        rows = inputdata.shape[1]
        columns = inputdata.shape[2]          
    else:          
        nints = 1
        ngroups = 1
        rows = inputdata.shape[0]
        columns = inputdata.shape[1]
        
    # Abort if the input data is not of the type expected
    if rows == 0 and columns == 0:
        # The input data does not contain a SCI data array
        strg = "Input file \'%s\' does not contain a SCI array. " % inputfile
        strg += "Are you sure this is level 1b exposure data?"
        raise TypeError(strg)
            
    # Add the reference output array, if it exists.
    if (input_model.refout is not None) and (len(input_model.refout) > 0):
        # Deduce the size of the reference output
        if input_model.refout.ndim >= 4:
            refrows = input_model.refout.shape[2]
            refcolumns = input_model.refout.shape[3]
            refaxis = 2
        elif input_model.refout.ndim >= 3:
            refrows = input_model.refout.shape[1]
            refcolumns = input_model.refout.shape[2]
            refaxis = 1       
        else:          
            refrows = input_model.refout.shape[0]
            refcolumns = input_model.refout.shape[1]
            refaxis =0

        new_refrows = int(refcolumns * refrows / columns)
        new_refcolumns = columns
        newshape = list(input_model.refout.shape)
        newshape[-2] = new_refrows
        newshape[-1] = new_refcolumns
        input_model.refout.shape = newshape
        newdata = np.concatenate( (inputdata, input_model.refout ), axis=refaxis )
        rows += new_refrows
    else:
        newdata = inputdata
        refrows = 0
        refcolumns = 0
        new_refrows = 0
        new_refcolumns = 0
                         
    readpatt = input_model.meta.exposure.readpatt
            
    output_model = ExposureData(rows, columns, ngroups, nints, readpatt)
    output_model.set_exposure(newdata)
    LOGGER.info("Data quality data is not copied during this conversion.")
            
    # Copy the metadata, keeping the same FITS keywords.
    # NOTE: No attempt is made to translate between different metadata standards.
    newheader = pyfits.Header()
    fits_dict = input_model.fits_metadata_dict()
    for kwd in (list(fits_dict.keys())):
        (fitshdu, fitskw, fitscomment) = fits_dict[kwd]
        if fitshdu == 'PRIMARY':
            try:
                newvalue = input_model.get_fits_keyword(fitskw, fitshdu)
                if newvalue is not None:
                    # DATE fields need to be explicitly converted to a string
                    if 'DATE' in fitskw:
                        newvalue = str(newvalue)
                    newheader[fitskw] = (newvalue, fitscomment)
            except (ValueError, TypeError) as e:
                strg = " Failed to convert metadata item with "
                strg += "FITS keyword \'%s\'. " % fitskw
                strg += "\n\t" + str(e)
                strg += " Value set to empty string."
                LOGGER.warning(strg)
                newheader[fitskw] = ('', '<-FAILED ' + fitscomment)
    # Define some explicit Level1b to DHAS header conversions
    newheader['NINT'] = nints
    newheader['NGROUP'] = ngroups
    try:
        newheader['READOUT'] = input_model.get_fits_keyword('READPATT')
    except KeyError:
        LOGGER.warning("READPATT keyword not found. READOUT=UNKNOWN.")
        newheader['READOUT'] = 'UNKNOWN'
    newheader['ROWSTART'] = 1
    newheader['ROWSTOP'] = rows - new_refrows
    if (input_model.refout is not None) and (len(input_model.refout) > 0):
        newheader['REFIMG'] = 'T'
        newheader['NREFIMG'] = new_refrows
    else:
        newheader['REFIMG'] = 'F'
    output_model.set_header(newheader)

    hstrg = "Created from level1b file %s" % inputfile
    output_model.add_fits_history(hstrg)
    hstrg = "Converted by SCASim convert_exposure_data script."
    output_model.add_fits_history(hstrg)
    output_model.save(outputfile, fileformat='FITSWriter', overwrite=overwrite)
    

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile"
    parser = optparse.OptionParser(usage)
    
#     # Optional arguments (long option strings only).
    parser.add_option("", "--convert", dest="convert", type="string",
                     default='level1b_to_dhas', help="Conversion choice"
                     )

    # Boolean options (short and long option strings).
    parser.add_option("-d", "--debug", dest="debug", action="store_true",
                      help="Debugging mode"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-s", "--silent", dest="silent", action="store_true",
                      help="Silent mode"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite existing file"
                     )

    (options, args) = parser.parse_args()
    
    # Compulsory arguments
    try:
        inputfile = args[0]
        outputfile = args[1]
    except IndexError:
        print( help_text )
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)
        
    # Optional arguments.
    convert_str = options.convert
        
    # Boolean flags.
    overwrite = options.overwrite
    verb = options.verb
    debug = options.debug
    silent = options.silent
    
    # Set the verbosity level according to the --verbose and --silent
    # options. (Note that --debug wins over --verbose and --silent wins
    # over all the other options if they are provided together.)
    verbose = 1
    if verb: verbose = 2
    if debug: verbose = 4
    if silent: verbose = 0

    # Ensure the output file is defined as .fits type
    if '.fits' not in outputfile:
        outputfile += '.fits'

    # Choose a conversion option
    if 'dhas_to_level' in convert_str.lower():
        if verbose > 0:
            print( "Converting DHAS file %s to level1b file %s." % \
                   (inputfile, outputfile) )
        convert_dhas_to_level1b( inputfile, outputfile, verbose=verbose,
                                 overwrite=overwrite )
        
    else: # level1b_to_dhas
        if verbose > 0:
            print( "Converting level1b file %s to DHAS file %s." % \
                   (inputfile, outputfile) )
        convert_level1b_to_dhas( inputfile, outputfile, verbose=verbose,
                                 overwrite=overwrite )

    print( "New exposure data saved to %s." % outputfile )
