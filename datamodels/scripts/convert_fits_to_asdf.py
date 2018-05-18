#!/usr/bin/env python
#
# :History:
# 
# 14 Mar 2017: Created.
# 30 Jun 2017: meta.reffile schema level removed to match changes in the
#              JWST build 7.1 data models release. meta.reffile.type also
#              changed to meta.reftype. TYPE keyword replaced by DATAMODL.
# 18 May 2018: Changed deprecated logger.warn() to logger.warning().
#
# @author: Steven Beard (UKATC)
#
"""

Script `convert_fits_to_asdf` reads a data model stored in a FITS-format
file and rewrites the model to an ASDF-format file. The data model must
be one of the recognised data types.

The following command arguments are defined by position::

    inputfile[0]
        The path+name of the FITS file to be read. Compulsory.
    outputfile[1]
        The path+name of the ASDF file to be written. Optional,
        but defaults to the same name as the input file.

The command also takes the following options::

    --datatype <type-string>
        The name of the data type to be used to read the product.
        If specified, this option overrides the TYPE keyword
        contained in the input file.
    --verbose or -v
        Generate more verbose output.

"""

from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import os, sys, time

# Python logging facility.
import logging
logging.basicConfig(level=logging.INFO)   # Choose ERROR, WARN, INFO or DEBUG 
LOGGER = logging.getLogger("convert_fits_to_asdf") # Get a default parent logger

# Import all the data models that might be contained in the file.
#from miri.datamodels.miri_measured_model import MiriMeasuredModel
import miri.datamodels        

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts a FITS file containing a MIRI data product into "
    usage += "an ASDF file."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--datatype", dest="datatype", type="string",
                     default=None, help="Data type to use (overriding TYPE)"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
#     parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
#                       help="Overwrite the copy of the file if it already exists"
#                      )

    (options, args) = parser.parse_args()
    if args and len(args) > 0:
        inputfile = args[0]
        inputbase, ext = os.path.splitext(inputfile)
        if isinstance(ext, bytes):
            ext = ext.decode(sys.getfilesystemencoding())
        if not ext:
            # Add a .fits extension if missing
            inputfile = inputfile + '.fits'
        elif ext != '.fits':
            # Otherwise the extension is supposed to be .fits.
            LOGGER.warning("Input file is supposed to be of \'.fits\' type, " + \
                           "not \'%s\'." % ext)
        if len(args) > 1:
            outputfile = args[1]
            if '.asdf' not in outputfile:
                outputfile = outputfile + ".asdf"     
        else:
            outputfile = inputbase + ".asdf"
    else:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    verb = options.verb
#     overwrite = options.overwrite

    if options.datatype:
        # Use the data type specified
        datatype = str(options.datatype)
        LOGGER.info("Forcing the data model to be opened with type \'%s\'" % datatype)
    else:
        datatype = ''

    # Display the data model using the class derived from the
    # data type.
    LOGGER.info("Opening data model from \'%s\'" % inputfile)
    with miri.datamodels.open( init=inputfile, astype=datatype ) as datamodel:
        if verb:
            if hasattr(datamodel.meta, 'reftype'):
                datatype = datamodel.meta.reftype
                strg = "The data model is of type \'%s\'" % str(datatype)
            else:
                strg = "The data model class name is \'%s\'" % \
                    datamodel.__class__.__name__
            LOGGER.info(strg)
                    
        LOGGER.info("Saving data model to \'%s\'." % outputfile)
        datamodel.save( outputfile )
        del datamodel
