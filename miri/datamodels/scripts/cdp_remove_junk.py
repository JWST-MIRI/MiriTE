#!/usr/bin/env python
#
# :History:
# 
# 03 Oct 2014: Created.
# 27 May 2015: Replaced pyfits with astropy.io.fits
# 12 Oct 2018: Resurrected to remove unwanted housekeeping HDUs and metadata
#              from MIRI CDP files.
#
# @author: Steven Beard (UKATC)
#

"""

Script `cdp_remove_junk` reads a FITS file saved by the JWST data models and
removes unwanted housekeeping HDUs, such as METADATA and ASDF. The script
also has the option to remove unwanted header keywords (such as the
GWAXTILT and GWAYTILT keywords not applicable to MIRI data).

The following command arguments are defined by position:

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_nojunk"
        appended.

The command also takes the following options:

    --hdus_to_remove or -H
        A comma-separated list of HDUs to be removed. If this option is not
        specified 'METADATA,ASDF' will be assumed.
    --keywords_to_remove or -K
        A command-separated list of keywords to be removed from the primary
        HDU. If this option is not specified, no keywords are removed.
    --verbose or -v:
        Generate more output.
    --overwrite or -c:
        Overwrite any existing FITS file.

"""

import optparse
import os, sys, time

import numpy as np
import astropy.io.fits as pyfits

def remove_junk( infile, outfile, hdu_list=None, keyword_list=None,
                 overwrite=False, verbose=False ):
    """

    Removes unwanted HDUs and primary metadata keywords from the given FITS
    file and writes a new file.
    
    """
    # This returned flag records whether a change has been made.
    changed = False
    
    # Attempt to read the FITS file
    with pyfits.open(infile) as hdulist:
        # Extract the primary header
        hdulist = pyfits.open(infile)
        header = hdulist[0].header
#         print "Primary header keywords:", header.keys()

        # Remove the unwanted primary header keywords.
        if keyword_list is not None and keyword_list:
            if verbose:
                print("Removing obsolete header keywords: %s" % str(keyword_list))
            for unwanted in keyword_list:
                if unwanted.strip() in header:
                    changed = True
                    del header[unwanted]
        
        # Remove the unwanted HDUs.
        if hdu_list is not None and hdu_list:
            for unwanted in hdu_list:
                kill = None
                for hi in range(0, len(hdulist)):
                    if hdulist[hi].name == unwanted.strip():
                        kill = hi
                if kill is not None:
                    if verbose:
                        print("Removing %s HDU." % unwanted)
                    changed = True
                    del hdulist[kill]
            
        # Write the result to a new file (ensuring the original name
        # is preserved in the header)
        if changed:
            try:
                hdulist[0].header['FILENAME'] = infile
            except (KeyError, IndexError, AttributeError):
                pass
            hdulist.writeto(outfile, overwrite=overwrite)

    # Close the hdulist and tidy up.
    try:
        hdulist.close()
        del hdulist
    except Exception:
        pass

    return changed

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Removes removes unwanted housekeeping HDUs and metadata created "
    usage += "by the JWST data models."
    parser = optparse.OptionParser(usage)
    parser.add_option("-H", "--hdus_to_remove", dest="hdus_to_remove", type="string",
                      default=None, help="List of HDUs to be removed"
                     )
    parser.add_option("-K", "--keywords_to_remove", dest="keywords_to_remove", type="string",
                      default=None, help="List of metadata keywords to be removed"
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
            outputfile = inputfile + "_nojunk.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    if options.hdus_to_remove is not None and options.hdus_to_remove:
        hdus_to_remove = str(options.hdus_to_remove)
        hdu_list = hdus_to_remove.split(",")
    else:
        hdu_list = ['METADATA', 'ASDF']

    if options.keywords_to_remove is not None and options.keywords_to_remove:
        keywords_to_remove = str(options.keywords_to_remove)
        keyword_list = keywords_to_remove.split(",")
    else:
        keyword_list = []

    verb = options.verb
    overwrite = options.overwrite

    # Convert the file
    print("Converting \'%s\'" %  inputfile)
    if remove_junk(inputfile, outputfile, hdu_list=hdu_list,
            keyword_list=keyword_list, overwrite=overwrite, verbose=verb):
        print("Saved to \'%s\'" % outputfile)
    else:
        print("No changes made. No output file saved.")
