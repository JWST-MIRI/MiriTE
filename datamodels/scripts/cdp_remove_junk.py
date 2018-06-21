#!/usr/bin/env python
#
# :History:
# 
# 03 Oct 2014: Created.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Steven Beard (UKATC)
#

"""

Script `cdp_remove_junk` reads a converted CDP-3 format file and removes
the junk that the conversion process leaves behind.

* FIELD_DEF and GROUP_DEF tables are removed if present.
* Obsolete metadata keywords (SUBXSIZE, SUBYSIZE, SUBXSTRT, SUBYSTRT)
  removed

The following command arguments are defined by position:

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out"
        appended.

The command also takes the following options:

    --overwrite or -c:
        Overwrite any existing FITS file.

"""



import optparse
import os, sys, time

import numpy as np
import astropy.io.fits as pyfits

def remove_junk( infile, outfile, overwrite=False ):
    """

    Removes the junk from the given file and writes a new file.
    
    """
    try:
        # Attempt to read the FITS file and extract the primary header.
        hdulist = pyfits.open(infile)
        header = hdulist[0].header
#         print "Primary header keywords:", header.keys()

        # Remove the obsolete header keywords.
        to_be_removed = ['SUBXSTRT', 'SUBXSIZE', 'SUBYSTRT', 'SUBYSIZE']
        print("Removing obsolete header keywords: %s" % str(to_be_removed))
        for hunted in to_be_removed:
            if hunted in header:
                del header[hunted]
        
        # Find and remove the FIELD_DEF and GROUP_DEF tables.
        to_be_removed = ['FIELD_DEF', 'GROUP_DEF']
        for hunted in to_be_removed:
            kill = None
            for hi in range(0, len(hdulist)):
                if hdulist[hi].name == hunted:
                    kill = hi
            if kill is not None:
                print("Removing %s HDU." % hunted)
                del hdulist[kill]
            
        # Write the result to a new file
        hdulist.writeto(outfile, overwrite=overwrite)

    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    return

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Removes the junk left over from the CDP-2 to CDP-3 "
    usage += "conversion process."
    parser = optparse.OptionParser(usage)
    parser.add_option("-c", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the FITS file if it already exists"
                     )

    (options, args) = parser.parse_args()

    try:
        inputfile = args[0]
        if len(args) > 1:
            outputfile = args[1]
        else:
            outputfile = inputfile + "_dejunked.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    overwrite = options.overwrite

    # Convert the file
    print("Converting \'%s\'" %  inputfile)
    remove_junk(inputfile, outputfile, overwrite=overwrite)
    print("Saved to \'%s\'" % outputfile)
