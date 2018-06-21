#!/usr/bin/env python
#
# :History:
# 
# 21 Jun 2016: Created.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
#
# @author: Steven Beard (UKATC)
#
"""

Script `find_me_another` takes a calibration data product, or a file
obtained from the JWST calibration reference system, and finds another
example from within the MIRI CDP repository.

A typical use of this script would be to replace a calibration file
extracted from the CRDS with a more recent version from the MIRI CDP
repository (before the new version has had a chance to be copied to the
CRDS).

Note that matching CDP files are extracted from the ftp repository and
are saved in the local CDP cache. A copy of this file is written to
the name specified in the outputfile parameter. A copy is written ONLY
if the input file is a recognised CDP and the file extracted from the
ftp repository has a different version number.

The following command arguments are defined by position::

    inputfile[0]
        The path+name of the file to be read as the example.
        
    outputfile[1]
        The path+name of the new file to be written.
        Optional. If not given, defaults to <inputfile>_new.fits.

The command also takes the following options::

    --datatype <type-string>
        The name of the data type to be used to read the product.
        If specified, this option overrides the TYPE keyword
        contained in the input file.

    --cdprelease
        The CDP release from which the new CDP file is to be imported.
        Defaults to the latest release.

    --cdpversion
        The CDP version from which the new CDP file is to be imported.
        Defaults to the latest version.

    --cdpsubversion
        The CDP subversion from which the new CDP file is to be imported.
        Defaults to the latest subversion.
        
    --verbose or -v
        Print the new model before saving it.

    --overwrite or -o
        Overwrite any existing FITS file.

"""



import optparse
import sys, time
import warnings

import astropy.io.fits as pyfits

from miri.datamodels.cdplib import get_cdp

def get_cdp_metadata( filename  ):
    """
    
    Helper function which extracts the CDP metadata from the FITS header
    of a file.
    
    """
    hdulist = None
    try:
        hdulist = pyfits.open( filename )

        if hdulist is not None:
            header = hdulist[0].header
            header_keys = list(header.keys())
            if 'REFTYPE' in header or 'REFTYPE' in header_keys:
                # There is a new data type keyword in the header.
                datatype = header['REFTYPE']
            elif 'TYPE' in header or 'TYPE' in header_keys:
                # There is an old data type keyword in the header.
                datatype = header['TYPE']
            else:
                datatype = ''  
            if 'DETECTOR' in header or 'DETECTOR' in header_keys:
                # There is a detector keyword in the header.
                detector = header['DETECTOR']
            else:
                detector = 'ANY'
            if 'READPATT' in header or 'READPATT' in header_keys:
                # There is a detector keyword in the header.
                readpatt = header['READPATT']
            else:
                readpatt = 'ANY'
            if 'SUBARRAY' in header or 'SUBARRAY' in header_keys:
                # There is a detector keyword in the header.
                subarray = header['SUBARRAY']
            else:
                subarray = 'ANY'
            if 'CHANNEL' in header or 'CHANNEL' in header_keys:
                # There is a filter keyword in the header.
                channel = header['CHANNEL']
            else:
                channel = 'ANY'
            if 'BAND' in header or 'BAND' in header_keys:
                # There is a filter keyword in the header.
                band = header['BAND']
            else:
                band = 'ANY'
            if 'FILTER' in header or 'FILTER' in header_keys:
                # There is a filter keyword in the header.
                mirifilter = header['FILTER']
            else:
                mirifilter = 'ANY'
            if 'VERSION' in header or 'VERSION' in header_keys:
                # There is a CDP version keyword in the header.
                version = header['VERSION']
            else:
                version = ''
            
    except Exception as e:
        strg = "Failed to open FITS file, \'%s\'\n" % filename
        strg += "  %s: %s" % (e.__class__.__name__, str(e))
        raise IOError(strg)
    finally:
        if hdulist is not None:
            hdulist.close()
            del hdulist
    return (datatype, detector, readpatt, subarray, channel, band,
            mirifilter, version)


if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile [outputfile]\n"
    usage += "Finds another (usually more recent) version of any "
    usage += "MIRI calibration data product."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--datatype", dest="datatype", type="string",
                     default=None, help="Data type to use (overriding TYPE)"
                     )
    parser.add_option("", "--cdprelease", dest="cdprelease", type="string",
                     default=None, help="CDP release to be searched"
                     )
    parser.add_option("", "--cdpversion", dest="cdpversion", type="string",
                     default=None, help="CDP version to be searched"
                     )
    parser.add_option("", "--cdpsubversion", dest="cdpsubversion", type="string",
                     default=None, help="CDP subversion to be searched"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the copy of the file if it already exists"
                     )

    (options, args) = parser.parse_args()

    try:
        inputfile = args[0]
        if len(args) > 1:
            outputfile = args[1]
        else:
            outputfile = inputfile + "_new.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    verb = options.verb
    overwrite = options.overwrite

    if options.datatype:
        # Use the data type specified
        datatype = str(options.datatype)
        print("Forcing the data model to be opened with type \'%s\'" % datatype)
    else:
        datatype = ''
        
    if options.cdprelease:
        cdprelease = options.cdprelease
    else:
        cdprelease = None
    if options.cdpversion:
        cdpversion = options.cdpversion
    else:
        cdpversion = None
    if options.cdpsubversion:
        cdpsubversion = options.cdpsubversion
    else:
        cdpsubversion = None

    # Obtain metadata from the example data model.
    (datatype, detector, readpatt, subarray, channel, band,
        mirifilter, version) = \
        get_cdp_metadata( inputfile )
    
    if datatype:
        # Attempt to find an alternative version of this data model
        strg = "Searching for a " + str(datatype) + " CDP"
        if detector != 'ANY':
            strg += ", DETECTOR=" + str(detector)
        if readpatt != 'ANY':
            strg += ", READPATT=" + str(readpatt)
        if subarray != 'ANY':
            strg += ", SUBARRAY=" + str(subarray)
        if band != 'ANY' or channel != 'ANY':
            strg += ", BAND=" + str(band) + ", CHANNEL=" + str(channel)
        if mirifilter != 'ANY':
            strg += ", FILTER=" + str(mirifilter)
        strg += "..."
        print(strg)
        newmodel = get_cdp(datatype, model='FM', detector=detector,
            readpatt=readpatt, channel=channel,
            band=band, mirifilter=mirifilter, subarray=subarray,
            integration=None,
            cdprelease=None, cdpversion=None, cdpsubversion=None,
            ftp_host=None, ftp_path=None, ftp_user='miri',
            ftp_passwd='', local_path=None, cdp_env_name='CDP_DIR',
            miri_env_name='MIRI_ENV', fail_message=True)
        
        # If this has worked, check the new model is different from the original
        # and, if so, save the new model to a new file.
        if newmodel is not None:
            if verb:
                print(newmodel)
            newversion = newmodel.meta.version
            if newversion != version:
                print("Saving new version of this CDP to %s." % outputfile)
                newmodel.save( outputfile, overwrite=overwrite )
            else:
                strg = "New CDP has exactly the same version number (%s)." % \
                    newversion
                strg += " Not saved."
                print(strg)
        else:
            print("Matching CDP could not be found.")

        del newmodel
    else:
        print("Input file does not look like a MIRI CDP. Data type unknown.")

