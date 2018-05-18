#!/usr/bin/env python
#
# :History:
# 
# 13 Mar 2017: Created.
# 30 Jun 2017: meta.reffile schema level removed to match changes in the
#              JWST build 7.1 data models release. meta.reffile.type also
#              changed to meta.reftype. TYPE keyword replaced by DATAMODL.
# 26 Jun 2017: Added cdp_ftp_path parameter and set to default if not specified.
# 17 May 2018: Python 3: Converted dictionary keys return into a list.
#
# @author: Steven Beard (UKATC)
#
"""

Script `cdp_get_doc` obtains the documentation associated with a given CDP
data product, copying one or more files to the local cache if necessary.

Typically, only one document should match a given CDP, but there might be
several issues of the same document or a set of documents on the same
subject. The optional "issue" parameter can be used to narrow the search to
a document whose filename contains a particular keyword or issue string
(e.g. "Issue1.4"). If left blank, all issues of a particular document are
returned, and the reader should check manually to see which one is more
relevant. (NOTE: It isn't possible to sort the issues automatically because,
unlike the CDPs themselves, each document uses a slightly different
versioning convention.)

The following command arguments are defined by position::

    inputfile[0]
        The path+name of the CDP file whose documentation is to be obtained.
        Compulsory.

The command also takes the following options::

    --username
        The username with which to access the CDP repository.
        The default is 'miri'
    --password
        The password with which to access the CDP repository.
        The default is blank, which will cause the script to
        prompt the user.
    --cdp_ftp_path
        A list of folders (or folders) on the SFTP host to be searched
        for CDP documents, consisting of a list of folder names separated by
        a ":" delimiter. Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'.
        Defaults to the standard CDP repository at Leuven.
    --issue
        An issue string which must be matched.
        The default is blank, which will match all issues.
    --datatype <type-string>
        The name of the data type to be used to read the CDP file.
        If specified, this option overrides the TYPE keyword
        contained in the input file.
    --force or -f
        Force an update even when the document already exists in the cache
    --verbose or -v
        Generate more verbose output.

"""

from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import sys, time

# Python logging facility.
import logging
logging.basicConfig(level=logging.INFO)   # Choose ERROR, WARN, INFO or DEBUG 
LOGGER = logging.getLogger("cdp_get_doc") # Get a default parent logger

# Import all the data models that might be contained in the file.
import miri.datamodels        
# Also import the MIRI CDP interface class
from miri.datamodels.cdplib import MiriCDPInterface

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile\n"
    usage += "Gets the document file(s) associated with "
    usage += "any MIRI calibration data product."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--username", dest="username", type="string",
                     default='miri', help="Username for CDP server"
                     )
    parser.add_option("", "--password", dest="password", type="string",
                     default=None, help="Password for CDP server"
                     )
    parser.add_option("", "--cdp_ftp_path", dest="cdp_ftp_path", type="string",
                     default='', help="Search path for imported CDPs"
                     )
    parser.add_option("", "--issue", dest="issue", type="string",
                     default=None, help="Document issue string to be matched"
                     )
    parser.add_option("", "--datatype", dest="datatype", type="string",
                     default=None, help="Data type to use (overriding TYPE)"
                     )
    parser.add_option("-f", "--force", dest="force", action="store_true",
                      help="Force an update from the repository"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )

    (options, args) = parser.parse_args()

    try:
        inputfile = args[0]
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)
        
    cdp_ftp_path = options.cdp_ftp_path
    if not cdp_ftp_path:
        cdp_ftp_path = MiriCDPInterface.FTP_PATH_DEFAULT

    force = options.force
    verb = options.verb
    if options.issue:
        issue = options.issue
    else:
        issue = ''
    if options.datatype:
        # Use the data type specified
        datatype = str(options.datatype)
        LOGGER.info("Forcing the data model to be opened with type \'%s\'" % datatype)
    else:
        datatype = ''

    miricdp = MiriCDPInterface(ftp_user=options.username,
                               ftp_passwd=options.password,
                               ftp_path=cdp_ftp_path,
                               logger=LOGGER)
    if verb:
        strg = "%s documents are available:\n" % len(miricdp.cdp_docs_available)
        LOGGER.info(strg)

    # Open the data model using the class derived from the data type.
    with miri.datamodels.open( init=inputfile, astype=datatype ) as datamodel:
        if verb:
            if hasattr(datamodel.meta, 'reftype'):
                datatype = datamodel.meta.reftype
                strg = "The data model is of type \'%s\'" % str(datatype)
            else:
                strg = "The data model class name is \'%s\'" % \
                    datamodel.__class__.__name__
            LOGGER.info(strg)

        # Get the history data associated with the data model,
        # search for the document string and extract the word
        # containing the document file name (assuming all file
        # names begin with 'MIRI').
        document_string = ''
        history_list = datamodel.get_history()
        for history in history_list:
            parts = list(history.keys())
            for key in parts:
                hist_string = str( history[key] )
                if 'DOCUMENT' in hist_string:
                    document_string = hist_string

        document_file = ''
        if document_string:
            words = document_string.split()
            for word in words:
                if 'MIRI' in word:
                    document_file = word

        if document_file:
            if verb:
                LOGGER.info("Document file declared in HISTORY is: " + document_file )
           
            cdp_docs = miricdp.get_cdp_doc( document_file, issue=issue,
                                            force_update=force )
            if len(cdp_docs) == 1:
               LOGGER.info("CDP document file \'%s\' copied to cache." % cdp_docs[0])
            elif len(cdp_docs) > 1:
               LOGGER.info("%d CDP document files copied to cache. Verify manually which one is the most relevant" % len(cdp_docs))
            else:
               LOGGER.error("No matching CDP documents found.")
        else:
           LOGGER.error("No HISTORY DOCUMENT record defined for this data model (or record not readable).")

        del datamodel
