#!/usr/bin/env python
#
# :History:
# 
# 05 Nov 2015: Created.
# 06 Nov 2015: Added option to change the USEAFTER and VERSION strings.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
# 30 Jun 2017: meta.reffile schema level removed to match changes in the
#              JWST build 7.1 data models release. meta.reffile.type also
#              changed to meta.reftype. TYPE keyword replaced by DATAMODL.
# 26 Apr 2018: Corrected exception raising syntax for Python 3.
#
# @author: Steven Beard (UKATC)
#
"""

Script `cdp_add_history` adds HISTORY records to an existing MIRI
calibration data product. The script can be used to add HISTORY records
to CDPs whose only problem is a missing HISTORY. The script can also
correct the USEAFTER and VERSION metadata.

The following command arguments are defined by position::

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out" appended.

The command also takes the following options::

    --history <history-string>
        A generic history string to be added.
        Defaults to 'Description of Reference File Creation'.
    --document <doc-name-string>
        Name of document describing the strategy and algorithms used to
        create file. Defaults to null string.
    --software <sw-descr-string>
        Description, version number, location of software used to create
        file. Defaults to null string.
    --dataused <data-used-string>
        Data used to create file. Defaults to null string.
    --differences <differences-string>:
        How is this version different from the one that it replaces?
        Defaults to null string.
    --useafter <use-after-date-string>
        Update the "USEAFTER" metadata to a new date, 'YYYY-MM-DD', where
        YYYY, MM and DD are the year, month and day of the month expressed
        numerically.
        If 'TODAY' is specified, today's date is used.
        If not specified, the USEAFTER metadata is not changed.
    --version <version-string>
        Update the "VERSION" metadata to a new string.
        If not specified, the VERSION metadata is not changed.
    --verbose or -v
        Generate more verbose output.
    --overwrite or -o
        Overwrite any existing FITS file.

"""

from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import sys, time

# Import all the data models that might be contained in the file.
#from miri.datamodels.miri_measured_model import MiriMeasuredModel
import miri.datamodels        

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Adds HISTORY records to "
    usage += "any MIRI calibration data product."
    parser = optparse.OptionParser(usage)
    defhist = 'Description of Reference File Creation'
    parser.add_option("", "--history", dest="history", type="string",
                     default=defhist, help="General history string"
                     )
    parser.add_option("", "--document", dest="document", type="string",
                     default="", help="Name of document describing CDP"
                     )
    parser.add_option("", "--software", dest="software", type="string",
                     default="", help="Software used to create CDP"
                     )
    parser.add_option("", "--dataused", dest="dataused", type="string",
                     default="", help="Data used to create CDP"
                     )
    parser.add_option("", "--differences", dest="differences", type="string",
                     default="", help="How is this CDP different from previous?"
                     )
    parser.add_option("", "--useafter", dest="useafter", type="string",
                     default="", help="New USEAFTER date (YYYY-MM-DD), if needed"
                     )
    parser.add_option("", "--version", dest="version", type="string",
                     default="", help="New VERSION string, if needed"
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
            outputfile = inputfile + "_out.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    history = options.history
    document = options.document
    software = options.software
    dataused = options.dataused
    differences = options.differences
    useafter = options.useafter
    version = options.version
    verb = options.verb
    overwrite = options.overwrite

    # Open the data model using the class derived from the data type.
    with miri.datamodels.open( init=inputfile ) as datamodel:
        
        # Add the history records.
        if history:
            datamodel.add_history(history)
        if document:
            hstrg = "DOCUMENT: " + document
            datamodel.add_history(hstrg)
        if software:
            hstrg = "SOFTWARE: " + software
            datamodel.add_history(hstrg)
        if dataused:
            hstrg = "DATA USED: " + dataused
            datamodel.add_history(hstrg)
        if differences:
            hstrg = "DIFFERENCES: " + differences
            datamodel.add_history(hstrg)
        if useafter:
            import time
            if hasattr(datamodel, 'meta'):
                if useafter == 'TODAY':
                    # Use today's date.
                    useafter = time.strftime('%Y-%m-%d')
                    datamodel.meta.useafter = useafter
                else:
                    # Check the string given matches the required format.
                    try:
                        test = time.strptime(useafter, '%Y-%m-%d')
                        datamodel.meta.useafter = useafter
                    except ValueError as e:
                        strg = "Invalid USEAFTER string, \'%s\'. " % useafter
                        strg += "Please use the date format \'YYYY-MM-DD\', "
                        strg += "e.g. \'2015-11-20\'."
                        raise ValueError(strg)
            else:
                strg = "Reference file metadata attributes missing from data model %s" % \
                    datamodel.__class__.__name_
                raise TypeError(strg)
        if version:
            datamodel.meta.version = version
            
        if verb:
            print(datamodel)
            print(datamodel.get_history_str())

        datamodel.save( outputfile, overwrite=overwrite)
        print("Data saved to %s\n" % outputfile)
            
        del datamodel
