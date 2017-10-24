# -*- coding: utf-8 -*-

"""

This module contains utilities for retrieving calibration data products
from their SFTP site. CDPs are stored in a local cache so they only need
to be retrieved once.

:Reference:

http://miri.ster.kuleuven.be/bin/view/Internal/CalDataProducts

:History:

25 Oct 2013: Created
05 Nov 2013: First version released
14 Feb 2014: Correct a problem where get_cdp() could download the entire
             cache if a matching file could not be found.
17 Feb 2014: Made some modifications to the regular expression string to
             help find CDPs.
20 Feb 2014: Improved the reliability of CDP searches by using global
             substring matching instead of regular expressions (which can
             cope with some of the CDPs being named in an inconsistent order).
             A regular expression search is still used to filter the
             version number. Also corrected a problem in the regular
             expression search - subarray names and filters can contain
             numbers.
21 Feb 2014: Modified to ensure a MiriCDPInterface object can be created
             even without a FTP connection. In this case, the class falls
             back to looking at the local cache only. Username and password
             options added for FTP access.
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
             Both sets of names are included here for backwards
             compatibility.
25 Sep 2014: Added extra channel and band settings.
02 Oct 2014: Global constants moved to util.py to prevent circular
             references.
08 Sep 2015: Made compatible with Python 3.
09 Sep 2015: Brought up to date with the restructuring of the MIRI CDP
             repository and the file naming convention for CDP-4.
             The CDP retrieval functions are working again.
25 Sep 2015: MiriCDPInterface class made a singleton, so the ftp connection
             is attempted only once per session. The local cache is now located
             by checking for the CDP_DIR environment variable before
             referencing MIRI_ENV. Corrected bug where local_path was not
             used when specified explicitly.
07 Oct 2015: Made exception catching Python 3 compatible.
07 Dec 2015: Allow the caller to specify a different logger.
23 Mar 2016: Only issue a warning message when there are multiple CDP
             candidates which differ in more than just the version number.
27 Apr 2016: Prompt for the ftp password when the link is first opened.
28 Apr 2016: Use getpass module to prompt for password without echoing.
04 May 2016: Mask the password string before logging.
05 May 2016: Added flat-field availability check. getpass does not work
             properly under Windows - raw_input used instead.
06 May 2016: Make get_cdp failure message optional.
06 Jun 2016: Added cdp_version_decode function. Corrected a bug where
             miri_env_name ended up None.
11 Jul 2016: Workaround for "IOError: [Errno ftp error] 200 TYPE is now 
             8-bit binary" encountered when using urllib on Linux: 
             explicitly reloading the urllib module prior to each file 
             retrieval.
08 Sep 2016: Do not clear the ftp password after a connection failure
             when that password has been given explicitly.
             Do not match the P750L filter when requesting a CDP for
             ANY filter for the MIRIMAGE detector.
             Added tests for simulator CDP files.
28 Oct 2016: Corrected a regular expression search which didn't expect
             a CDP type string to contain a "-" character.
01 Oct 2016: Added an extra filter which ensures that band strings don't
             accidentally match with the detector names (e.g. LONG matching
             with MIRIFULONG rather than with 34LONG).
11 Jan 2017: Include the mirifilter='GENERIC' option, which can be used to
             avoid CDPs designed for one specific filter. Corrected potential
             issue where function parameters were not passed by keyword.
             Changed the search logic so subarray=None is the equivalent of
             subarray='FULL' and not subarray='ANY'.
19 Jan 2017: Replaced use of ftp and urllib with pysftp. Default username
             for simulator CDPs changed from miriuser to cdpuser.
             Encapsulate the sftp connection in _open and _close methods.
             Do not count files within sub-directories of the local cache.
24 Jan 2017: The list of CDP filenames is now sorted so that pre-release
             versions, such as 5B and 6B come before full release versions
             such as 05 and 06. Regular expression match corrected so that
             CDP release numbers have 1 to 2 digits, not 2 digits.
23 Feb 2017: Check whether pysftp supports CnOpts before using it.
09 Mar 2017: Add documentation on setting ftp_host to 'LOCAL'.
13 Mar 2017: Look for documentation as well as calibration data.
17 Mar 2017: Corrected some documentation typos.
20 Apr 2017: Added a timeout option. Make hostkey checking the default
             if it is available and give a warning if a hostkey is
             not found.
09 May 2017: Major change to allow a search path of FTP folders to be
             provided. CDPFolder class added to manage each CDP folder.
             Tests updated.
07 Jul 2017: Give a warning if the CDP cache contains an empty file.
21 Jul 2017: Configure simulator tests to use full FTP search path.
26 Jul 2017: Top-level API changed so that MiriCDPInterface search functions
             return a list of filenames and a corresponding list of folders.
             Explicitly convert string parameters to string.
04 Oct 2017: Recent updates broke the version number sorting update made
             on 24 Jan 2017. Now corrected.
12 Oct 2017: Make the detector parameter mandatory in get_cdp(), since it
             is a primary key into the CDP_DICT dictionary.
24 Oct 2017: Corrected persistence problem where ftp-host remained 'LOCAL'
             even when changed back. 

Steven Beard (UKATC), Vincent Geers (UKATC)

"""

# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

from astropy.extern import six
import os
import re
import time
import sys, getpass
import copy

# Python utilities for accessing the Ftp repository.
import pysftp
from paramiko import SSHException

# Python logging facility.
import logging
logging.basicConfig(level=logging.INFO)   # Turn off verbose paramiko messages
LOGGER = logging.getLogger("miri.cdplib") # Get a default parent logger
# Logging level for the CDP classes
LOGGING_LEVEL = logging.INFO # Choose ERROR, WARN, INFO or DEBUG

# Import CDP utility functions and CDP dictionary.
from miri.datamodels.util import MIRI_MODELS, MIRI_DETECTORS, \
    MIRI_SETTINGS, MIRI_READPATTS, MIRI_SUBARRAYS, MIRI_CHANNELS, \
    MIRI_BANDS, MIRI_FILTERS, get_data_class
from miri.datamodels.cdp import CDP_DICT

# List all public classes and global functions here.
__all__ = ['get_cdp', 'MiriCDPInterface']

#
# (1) Global functions
#
def _criteria_string(cdptype, model='FM', detector=None,
                     readpatt=None, channel=None, band=None,
                     mirifilter=None, subarray='FULL', integration=None,
                     cdprelease=None, cdpversion=None, cdpsubversion=None):
    """
    
    Helper function which packages the CDP search criteria into a readable string
    
    """
    strg = cdptype
    if detector and (detector != 'ANY'):
        strg += ", detector=" + detector
    if readpatt and (readpatt != 'ANY'):
        strg += ", readpatt=" + readpatt
    if ((channel and (channel != 'ANY')) or \
        (band and (band != 'ANY'))):
        strg += ", channel/band="
        if channel and (channel != 'ANY'):
            strg += str(channel)
        if band and (band != 'ANY'):
            strg += str(band)
    if mirifilter and (mirifilter != 'ANY'):
        strg += ", filter=" + mirifilter
    if subarray and (subarray != 'ANY') and (subarray != 'FULL'):
        strg += ", subarray=" + subarray
    if integration:
        strg += ", integration=" + str(integration)
    if cdprelease or cdpversion or cdpsubversion:
        strg += ", release="
        if cdprelease:
            strg += str(cdprelease)
        if cdpversion or cdpsubversion:
            strg += "."
            if cdpversion:
                strg += str(cdpversion)
            if cdpsubversion:
                strg += "." + str(cdpsubversion)
    return strg

def _diff_count(string1, string2):
    """
    
    Count the number of characters by which two strings differ.
    
    """
    assert isinstance(string1, str)
    assert isinstance(string2, str)
    if string1 == string2:
        return 0
    minlen = min(len(string1), len(string2))
    diffcount = abs(len(string1) - len(string2))
    for ii in range(0,minlen):
        if string1[ii] != string2[ii]:
            diffcount += 1
    return diffcount

def _diff_count_list( string_list ):
    """
    
    Return the average character difference between pairs of strings in a list
    
    """
    assert isinstance(string_list, (tuple,list))
    if len(string_list) < 2:
        return 0
    diffsum = 0
    previous = None
    for entry in string_list:
        if previous is not None:
            diffsum += _diff_count(previous, entry)
            previous = entry
    return diffsum // (len(string_list) - 1)

def cdp_version_decode( cdp_version ):
    """
    
    Decode a CDP version string of the form 'xx.yy.zz' into
    xx -> cdprelease, yy --> cdpversion and zz -> cdpsubversion.
    Strings of the form 'xx.yy' or 'xx' are also decoded (with the
    missing version numbers returned as None).
    
    If an empty or invalid string is provided, all three version
    numbers are returned as None.

    :Parameters:
    
    cdp_version: str
        CDP version string, of the form 'xx.yy.zz', 'xx.yy' or 'xx'.

    :Returns:
    
    xxyyzz: tuple of 3 str
        (xx, yy, zz)
    
    """
    # Reject a null version code, a non-string or an empty string.
    if cdp_version is None or \
       not isinstance(cdp_version, (str,unicode)) or \
       len(cdp_version) < 1:
        return (None, None, None)
    
    # Split the version code into words separated by '.'.
    version = cdp_version.split('.')
    nwords = len(version)
    cdprelease = version[0]
    if nwords > 1:
        cdpversion = version[1]
        if nwords > 2:
            cdpsubversion = version[2]
        else:
            cdpsubversion = None   
    else:
        cdpversion = None
        cdpsubversion = None
    
    return (cdprelease, cdpversion, cdpsubversion)

#+++ MAIN FUNCTION +++
def get_cdp(cdptype, detector, model='FM', readpatt='ANY', channel='ANY',
            band='ANY', mirifilter='ANY', subarray='FULL', integration=None,
            cdprelease=None, cdpversion=None, cdpsubversion=None,
            ftp_host=None, ftp_path=None, ftp_user='miri', ftp_passwd='',
            timeout=None, local_path=None, cdp_env_name='CDP_DIR',
            miri_env_name='MIRI_ENV', logger=LOGGER, fail_message=True):
    """
    
    Get a calibration data product matching the specified criteria.
    The minimum criterion is the type of product required. If other
    criteria are not specified they are assumed to be wild cards
    which match everything. If an explicit release and version number
    is not given, the most recent version will be obtained.

    NOTE: This function returns ONE data object, even if more than one
    CDP file matches the criteria given. The last entry in an
    alphanumerically sorted list will be returned. It is recommended that
    all the required parameters are specified.

    If no CDP file could be matched, the function returns None.
    
    :Parameters:
    
    cdptype: str
        Data type of CDP to be obtained. Must be one of the defined
        MIRI CDP types. For example, 'MASK', 'DARK', 'LASTFRAME',
        'DISTORTION', 'DROOP', 'PIXELFLAT', 'FRINGE', 'PHOTOM', 'IPC',
        'COLCORR', 'FLUX', 'JUMP', 'LATENT', 'LINEARITY', 'SATURATION',
        'PSF', 'STRAY', 'TRACORR', 'WAVCORR', 'TELEM', etc....
    detector: str
        The MIRI detector required ('MIRIMAGE', 'MIRIFUSHORT',
        'MIRIFULONG' or 'ANY').
        The old names ('IM', 'SW' or 'LW') may be used to find CDPs
        prior to the CDP-3 release.
    model: str (optional)
        The MIRI model required ('VM', 'FM' or 'JPL').
        Defaults to 'FM'.
    readpatt: str (optional)
        The MIRI readout pattern required ('FAST', 'SLOW', or 'ANY').
        Defaults to 'ANY', which will match any readout pattern.  
    channel: str (optional)
        The MIRI MRS channel required ('1', '2', '3', '4', '12',
        '34' or 'ANY').
        Valid only for MRS data, when detector is 'MIRIFUSHORT' or
        'MIRIFULONG'.
        Defaults to 'ANY', which will match any channel.  
    band: str (optional)
        The MIRI MRS band required ('SHORT', 'MEDIUM', 'LONG',
        'SHORT-MEDIUM', 'SHORT-LONG', 'MEDIUM-SHORT', etc..., or 'ANY').
        Valid only for MRS data, when detector is 'MIRIFUSHORT' or
        'MIRIFULONG'.
        Defaults to 'ANY', which will match any band.  
    mirifilter: str (optional)
        The MIRI filter required ('F560W', 'F770W', 'F1000W', 'F1130W',
        'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W', 'F2550WR',
        'F1065C', 'F1140C', 'F1550C', 'F2300C', 'P750L', 'FLENS', 'FND',
        'OPAQUE', 'ANY' or 'GENERIC').
        Valid only for imager data, when detector is 'IM'.
        Defaults to 'ANY', which will match any filter.
        Explicitly specifying 'GENERIC' will only match CDPs which do not
        specify any filter and therefore work with any filter (as opposed
        to 'ANY', which can match a CDP which is specific to one filter only).
    subarray: str (optional)
        The MIRI subarray required ('FULL', 'MASK1140', 'MASK1550',
        'MASK1065', 'MASKLYOT', 'BRIGHTSKY', 'SUB256', 'SUB128',
        'SUB64' or 'SLITLESSPRISM').
        Defaults to 'FULL', which will match full frame data.
    integration: int (optional)
        When a CDP is split between separate files for each integration,
        this parameter can be used to specify a particular integration.
        Valid for 'DARK' and 'LINEARITY' data only.
        Defaults to any integration.
    cdprelease: int or str (optional)
        A string or integer describing the CDP release number required.
        If not specified, the latest release will be used.
    cdpversion: int or str (optional)
        A string or integer describing the CDP format version number required.
        If not specified, the latest version will be used.
    cdpsubversion: int or str (optional)
        A string or integer describing the CDP subversion number required.
        If not specified, the latest subversion will be used.

    :AdditionalParameters:

    ftp_host: str (optional)
        The address of the machine hosting the MIRI SFTP repository.
        Defaults to the MIRI CDP repository at Leuven.
        Setting this to 'LOCAL' will force MIRI CDPs to be obtained from
        the local cache.
    ftp_path: str (optional)
        The path to the folder (or folders) on the SFTP host where the
        MIRI CDPs are held. A search path consisting of more than one
        folder may be specified, separated by a ":" delimiter. These
        folders will be checked in the order given, and the first CDP
        matching the given criteria used.
        Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
        Defaults to the default CDP repository at Leuven.
    ftp_user: str (optional)
        A username with which to access the ftp site.
        Defaults to 'miri'.
    ftp_passwd: str (optional)
        A password with which to access the ftp site.
    timeout: float (optional)
        A new connection timeout in seconds.
        If not specified, the default SFTP connection timeout is used.
    local_path: str (optional)
        A local file path giving the MIRI CDP folder.
        If not given, the local path will be obtained from the
        environment variable specified in the cdp_env_name parameter,
        or failing that from the miri_env_name parameter.
        If local_path is not defined, and none of the named environment
        variables are defined, the local path defaults to '.'.
        Examples: '.',  '/home/me/MIRI_DHAS/MPipeline/Cal".
    cdp_env_name: str (optional)
        The name of the environment variable defining the MIRI CDP path.
        Defaults to 'CDP_DIR'.
        If the named environment variable is not defined, the miri_env_name
        parameter is used to obtain the top-level MIRI environment name.
    miri_env_name: str (optional)
        The name of the environment variable containing the top level MIRI
        path. CDPs are assumed to be contained in a '/CDP' subdirectory.
        Defaults to 'MIRI_ENV'.
    logger: Logger object (optional)
        A Python logger to handle the I/O. This parameter can be used
        by a caller to direct the output to a different logger, if
        the default defined by this module is not suitable.
    fail_message: boolean (optional)
        When True (the default) a failure message is logged when a
        CDP cannot be found and None is returned.
        Set to False when using get_cdp to try different options and
        you don't want a message logged for each try.
        
    :Environment:
  
    cdp_env_name (CDP_DIR):
        Environment variable defining the MIRI CDP path.
        See the cdp_env_name parameter.
    miri_env_name (MIRI_ENV)
        Environment variable defining the top level MIRI path.
        See the miri_env_name parameter.

    :DataStorage:
    
    When a CDP is obtained using this function, a local copy is stored in
    a cache directory defined by the local_path, cdp_env_name and miri_env_name
    parameters. local_path overrides cdp_env_name, which overrides
    miri_env_name. See the MiriCDPInterface documentation.
    
    :Returns:
    
    dataproduct: JWST data product
        A data product containing the specified CDP. Returns None
        if the data product could not be found or accessed.
    
    """
    # Get a logging object
    mylogger = logger.getChild("get_cdp")
    
    # The cdptype must be one of the available CDP types
    if not cdptype in CDP_DICT:
        # Not a recognised CDP type.
        strg = "Data type \'%s\' is not a recognised MIRI CDP.\n" % cdptype
        strg += "It must be one of: "
        start = True
        for key in CDP_DICT.keys():
            if start:
                start = False
            else:
                strg += ", "
            strg += "\'%s\'" % key
        raise TypeError(strg)
    
    # Access the MIRI CDP repository through a CDP interface object.
    # NOTE: The MiriCDPInterface is a singleton, so the class is created
    # once, per session, even if get_cdp is called many times.
    CDPInterface = MiriCDPInterface(ftp_host=ftp_host, ftp_path=ftp_path,
                                    ftp_user=ftp_user, ftp_passwd=ftp_passwd,
                                    timeout=timeout, local_path=local_path,
                                    cdp_env_name=cdp_env_name,
                                    miri_env_name=miri_env_name, logger=logger)
    
    # Refresh the interface if any parameters are different from when the
    # class was first created (necessary when the class is a singleton).
    CDPInterface.refresh(ftp_host=ftp_host, ftp_path=ftp_path,
                         ftp_user=ftp_user, ftp_passwd=ftp_passwd,
                         timeout=timeout, local_path=local_path,
                         cdp_env_name=cdp_env_name,
                         miri_env_name=miri_env_name)
    
    # Get the name of a CDP file matching the specified criteria.
    (filename, ftp_path) = CDPInterface.match_cdp_latest(cdptype, model=model,
                    detector=detector, readpatt=readpatt, channel=channel,
                    band=band, mirifilter=mirifilter, subarray=subarray,
                    integration=integration, cdprelease=cdprelease,
                    cdpversion=cdpversion, cdpsubversion=cdpsubversion)         
    
    if filename:
        # Update the local CDP cache to make sure it contains the specified file,
        # and obtain the local file path and name.
        local_filename = CDPInterface.update_cache(filename, ftp_path)

        # Read the contents of the file into a new data model, using the class
        # derived associated with the data type, detector and filter.
        kwlist = [cdptype]
        if not detector or detector == 'ANY':
            # Ambiguous detector. Try to figure it out from the other parameters.
            if mirifilter:
                detector = 'MIRIIMAGE'
            elif miriband:
                if '1' in miriband or '2' in miriband:
                    detector = 'MIRIFUSHORT'
                elif '3' in miriband or '4' in miriband:
                    detector = 'MIRIFULONG'
            strg = "Detector specification is ambiguous. Trying detector=%s" % \
                str(detector)
            mylogger.warn(strg)

        if detector:
            kwlist.append(detector)
        if mirifilter:
            kwlist.append(mirifilter)
        data_class = get_data_class(kwlist)
        strg = "Reading \'%s\' model from \'%s\'" % (cdptype, local_filename)
        mylogger.info(strg)
        datamodel = data_class( init=local_filename )
    else:
        if fail_message:
            criteria = _criteria_string(cdptype, model=model,
                detector=detector, readpatt=readpatt, channel=channel,
                band=band, mirifilter=mirifilter, subarray=subarray,
                integration=integration, cdprelease=cdprelease,
                cdpversion=cdpversion, cdpsubversion=cdpsubversion)
            strg = "The criteria given (%s) did not match any CDP files." % criteria
            mylogger.error(strg)
        datamodel = None

    # If successful, return the new data model.
    return datamodel

#
# (2) Classes
#
# A metaclass which only allows one instance of a class to exist at a time.
class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = \
                super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class MiriCDPFolder(object):
    """
    
    A class which manages a FTP folder containing MIRI Calibration Data
    Products. It contains search functions for locating particular CDPs
    within that folder.
    
    The class is designed to be contained within a MiriCDPInterface object.
    Objects exists for the lifetime of a particular STFP/folder combination.
    If the SFTP connection is closed, or if the FTP path changes, the
    object is deleted. 
    
    :Parameters:
    
    sftp: pysftp Connection object
        An open SFTP connection, managed by pysftp.
        None if local host used.
    ftp_path: str
        A file path on the FTP repository.
        Only used if an STFP connection is open.
    cdp_dir: str
        A path to the local CDP folder.
                
    """
                                             
    def __init__(self, sftp, ftp_path, cdp_dir, logger=LOGGER):
        """
        
        Initialises the MiriCDPFolder class.
        
        Parameters: See class doc string.

        """
        # Get a Python logger object
        self.logger = logger.getChild(self.__class__.__name__)
        self.logger.setLevel(level=LOGGING_LEVEL)
        self.logger.debug("Class %s created for ftp_path=\'%s\', cdp_dir=\'%s\'."  % \
                          (self.__class__.__name__, str(ftp_path), str(cdp_dir)))

        self.sftp = sftp
        self.ftp_path = str(ftp_path)
        self.cdp_dir = str(cdp_dir)
        
        # Initialise the list of available CDPs
        self.cdp_files_available = []
        self.cdp_docs_available = []
        self.update_cdp_list()

    def update_cdp_list(self):
        """
        
        Helper function to update the current list of CDPs.
        The function connects to the MIRI CDP SFTP repository and
        requests a file listing, which is then converted to a list.
        Only files whose names begin with MIRI_ are included.
        
        If a connection cannot be made to the SFTP repository,
        the function looks at the files contained in the local
        cache. NOTE: In this mode the class is only made aware
        of files that have already been downloaded to the cache.
        If there are frequent SFTP problems, manually copy all
        the available CDPs to the local cache before using this
        utility.

        :Parameters:
    
        timeout: float (optional)
            If provided, a new connection timeout in seconds.
            By default, the connection timeout is not changed.
        
        """
        if self.cdp_files_available:
            del self.cdp_files_available
        if self.cdp_docs_available:
            del self.cdp_docs_available
        ftp_list = []
        msg = ''
        
        if self.sftp is not None:
            try:
                self.sftp.chdir(self.ftp_path)
            except IOError, e:
                strg = "IOError: Failed to change directory to FTP folder \'%s\'\n" % self.ftp_path
                strg += "  %s" % str(e)
                raise IOError(strg)
            ftp_list = self.sftp.listdir()
            self.sftp.chdir('/')
            self.logger.debug("CDP list updated for ftp_path=\'%s\'." % self.ftp_path)
        else:
            strg = "No ftp host specified. Using local cache.\n  "
            self.logger.warning(strg)

            # Make sure the local cache exists.
            abspath = os.path.abspath(self.cdp_dir)
            if os.path.isdir(abspath):
                # Walk through the local cache and find file names.
                for (dirpath, dirnames, filenames) in os.walk(abspath):
                    for filename in filenames:
                        ftp_list.append(filename)
                    # Break to ensure only the top-level directory is included.
                    break
            else:
                strg = "Local cache \'%s\' does not exist.\n" % abspath
                strg += "  If you can't access the SFTP site, please create "
                strg += "cache directory and copy CDP files manually."
                raise IOError(strg)

        # Only add names matching MIRI CDPs (which start 'MIRI_') to the list.
        # Restrict the list to files containing '.fits' to exclude the PDF
        # and text documentation.
        self.cdp_files_available = []
        for ftp_file in ftp_list:
#             self.logger.debug("Testing ftp_file=%s", ftp_file)
            if re.match("^MIRI_", ftp_file) is not None and \
               '.fits' in ftp_file:
                self.cdp_files_available.append(ftp_file)

        # Give a warning if no files have been found. This will cause
        # all subsequent attempts to match or get a CDP file to fail.
        if not self.cdp_files_available:
            self.logger.warn("No CDP files available!")

        # Find the documents associated with the calibration files
        self.cdp_docs_available = []
        for ftp_file in ftp_list:
#             self.logger.debug("Testing ftp_file=%s", ftp_file)
            for doc_type in ('.pdf', '.doc', '.txt'):
                if re.match("^MIRI", ftp_file) is not None and \
                   doc_type in ftp_file:
                    self.cdp_docs_available.append(ftp_file)
                        
    def _filter_regexp(self, input_list, match_string, flags=0):
        """
        
        Filter a list of strings and return a shorter list containing
        only those matching a given regular expression.
        
        :Parameters:
        
        input_list: list of str
            A list of strings to be filtered.
        match_string: raw string
            A regular expression to be matched.
            For example, r'[A-Za-z_]*Bad[A-Za-z_]*.+fits' matches
            all bad pixel mask file names.
        flags: int (optional)
            Flags used to modify the regular expression search.
            For example flags=re.IGNORECASE will treat lowercase and
            uppercase matches equally.
            The default matches the exact case.
            
        :Returns:
        
        matched_list: list of str
            A list of strings matching the regular expression
        
        """
        matched_list = []
        for test_string in input_list:
            if re.match(match_string, test_string, flags=flags) is not None:
                matched_list.append(test_string)                
        return matched_list
                
    def match_cdp_regexp(self, match_string, flags=re.IGNORECASE):
        """
        
        Return a sorted list of all the CDP files contained in the FTP
        folder matching a given regular expression.
        
        :Parameters:
        
        match_string: raw string
            A regular expression to be matched.
            For example, r'[A-Za-z_]*Bad[A-Za-z_]*.+fits' matches
            all bad pixel mask file names.
        flags: int (optional)
            Flags used to modify the regular expression search.
            For example, the default flags=re.IGNORECASE will treat 
            lowercase and uppercase matches equally.
            Specific flags=0 to match the exact case.
            
        :Returns:
        
        matched_files: list of str
            A sorted list of matching filenames
        
        """
        self.logger.debug("Matching regexp: " + match_string)
        matched_files = self._filter_regexp(self.cdp_files_available,
                                             match_string, flags=flags)
                
        # Return a sorted list
        matched_files.sort()
        return matched_files

    def match_cdp_substrings(self, mustcontain=[], mustnotcontain=[]):
        """
        
        Return a sorted list of all the CDP files contained in the FTP
        folder whose names contain all the compulsory substrings and
        do not contain any forbidden substrings.
        
        :Parameters:
        
        mustcontain: list of str
            A list of substrings which must be matched.
        mustnotcontain: list of str
            A list of substrings which must NOT be matched.
            
        :Returns:
        
        matched_files: list of str
            A sorted list of matching filenames
        
        """
        self.logger.debug("Must contain: " + str(mustcontain))
        self.logger.debug("Must not contain: " + str(mustnotcontain))
        matched_files = []
        for cdp_file in self.cdp_files_available:
            cdp_file_u = cdp_file.upper()
            strings_matched = True
            for match_string in mustcontain:
                if match_string.upper() not in cdp_file_u:
                    strings_matched = False
                    break
            strings_avoided = True
            for match_string in mustnotcontain:
                if match_string.upper() in cdp_file_u:
                    strings_avoided = False
                    break
            if strings_matched and strings_avoided:
                matched_files.append(cdp_file)
                
        # Return a sorted list
        matched_files.sort()
        return matched_files

    def _encode_version(self, code):
        """
    
        Helper function to convert a specified code into a formatted
        version number. For example, 2 is formatted into '02'.
         
        :Parameters:
        
        code: int or str
            Version number to be formatted
            
        :Returns:
        
        formatted_version: str
            A formatted version number.
    
        """
        if code is None:
            # The version number is a place holder.
            version = '[0-9]{2}'
        elif isinstance(code, (int,float)):
            # Convert from integer to string, assuming 2 digits.
            version = '%.2d' % int(code)
        else:
            # Assume already a formatted string
            version = str(code)
        return version

    def match_cdp_filename(self, cdptype, model='FM', detector=None,
                           readpatt=None, channel=None, band=None,
                           mirifilter=None, subarray='FULL', integration=None,
                           cdprelease=None, cdpversion=None, cdpsubversion=None):
        """
        
        Find a list of available CDPs contained in the FTP folder matching the
        required criteria. The CDP file naming convention is assumed to be:
        
        MIRI_<model>_<detector>_<detsetng>_<readpatt>_<channelband_or_filter>_<subarray>_<reftype>_<version>.fits
        
        :Parameters:
    
        cdptype: str
            Data type of CDP to be obtained. Must be one of the defined
            MIRI CDP types. For example, 'MASK', 'DARK', 'LASTFRAME',
            'DISTORTION', 'DROOP', 'PIXELFLAT', 'FRINGE', 'PHOTOM', 'IPC',
            'COLCORR', 'FLUX', 'JUMP', 'LATENT', 'LINEARITY', 'SATURATION',
            'PSF', 'STRAY', 'TRACORR', 'WAVCORR', 'TELEM', etc....
        model: str (optional)
            The MIRI model required ('VM', 'FM' or 'JPL').
            Defaults to 'FM'.
        detector: str (optional)
            The MIRI detector required ('MIRIMAGE', 'MIRIFUSHORT',
            'MIRIFULONG' or 'ANY'). The old names ('IM', 'SW' or 'LW')
            may be used to find CDPs prior to the CDP-3 release.
            By default, matches all detectors.
        readpatt: str (optional)
            The MIRI readout pattern required ('FAST', 'SLOW' or 'ANY').
            By default, matches all readout patterns.  
        channel: str (optional)
            The MIRI MRS channel required ('1', '2', '3', '4', '12',
            '34' or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches all channels.  
        band: str (optional)
            The MIRI MRS band required ('SHORT', 'MEDIUM', 'LONG',
            'SHORT-MEDIUM', 'SHORT-LONG', 'MEDIUM-SHORT', etc..., or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches all bands.  
        mirifilter: str (optional)
            The MIRI filter required ('F560W', 'F770W', 'F1000W', 'F1130W',
            'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W', 'F2550WR',
            'F1065C', 'F1140C', 'F1550C', 'F2300C', 'P750L', 'FLENS', 'FND',
            'OPAQUE', 'ANY' or 'GENERIC').
            Valid only for imager data, when detector is 'IM'.
            By default, matches all filters.  
            Explicitly specifying 'GENERIC' will only match CDPs which do not
            specify any filter and therefore work with any filter (as opposed
            to 'ANY', which can match a CDP which is specific to one filter only).
        subarray: str (optional)
            The MIRI subarray required ('FULL', 'MASK1140', 'MASK1550',
            'MASK1065', 'MASKLYOT', 'BRIGHTSKY', 'SUB256', 'SUB128',
            'SUB64', 'SLITLESSPRISM').
            Defaults to 'FULL', which will match full frame data.
        integration: int (optional)
            When a CDP is split between separate files for each integration,
            this parameter can be used to specify a particular integration.
            Valid for 'Dark' data only.
            By default matches all integrations.
        cdprelease: int or str (optional)
            A string or integer describing the CDP release number required.
            If not specified, the latest release will be matched.
        cdpversion: int or str (optional)
            A string or integer describing the CDP format version number
            required.
            If not specified, the latest version will be matched.
        cdpsubversion: int or str (optional)
            A string or integer describing the CDP subversion number required.
            If not specified, the latest subversion will be matched.
        
        :Returns:
        
        filenames: list of str
            A sorted list of names of files matching the given criteria.
            Returns an empty list if no files are matched.
         
        """
        if self.sftp is None:
            self.logger.warn("Matching against local CDP cache only.")
        
        # Build a list of substrings that must be contained and must be avoided.
        match_strings = []
        avoid_strings = []
        
        # The filename always contains a model name.
        if model is not None and model in MIRI_MODELS:
            # Match the exact MIRI model
            match_strings.append(model)

        # The model name is always followed by a detector name, unless
        # the model name is 'JPL'
        if model != 'JPL':
            if detector is not None and detector in MIRI_DETECTORS:
                # Match the exact detector
                match_strings.append(detector)
        
        # The readout pattern is optional, and a CDP valid for any
        # pattern is specified by missing it out completely.
        if readpatt is not None and readpatt != 'ANY':
            match_strings.append(readpatt)
               
        # The filter name is optional, and a CDP valid for any
        # filter is specified by missing it out completely.
        # Specifying 'GENERIC' will avoid CDPs designed for specific filters.
        if mirifilter is not None and mirifilter != 'ANY' and \
           mirifilter != 'GENERIC':
            match_strings.append(mirifilter)
        elif mirifilter is not None and mirifilter == 'GENERIC':
            for filt in MIRI_FILTERS:
                avoid_strings.append(filt)
            
        # If an imager CDP is needed without specifying a filter,
        # explicitly exclude the LRS CDPs, for which either the
        # 'P750L' filter or the 'SLITLESSPRISM' subarray would
        # have been specified explicitly.
        if detector == 'MIRIMAGE' and mirifilter == 'ANY' and \
           not subarray == 'SLITLESSPRISM':
            avoid_strings.append('P750L')

        # The channel and band names are optional, and a CDP valid for any
        # channel or band is specified by missing either out completely.
        channel_band = ''
        if channel is not None and channel != 'ANY':
            channel_band += str(channel)
        if band is not None and band != 'ANY':
            channel_band += band
        if channel_band:
            match_strings.append(channel_band)

        # The subarray is optional. A full frame CDP is represented by missing
        # it out completely. If a specific subarray is specified, only that
        # subarray should be matched. If FULL frame data are needed, specific
        # subarrays should be explicitly avoided. Only the subarray='ANY'
        # option will match any subarray.
        if subarray is not None and subarray != 'FULL' and \
           subarray != 'ANY' and subarray != 'GENERIC':
            match_strings.append(subarray)
        if subarray == None or subarray == 'FULL' or subarray == 'GENERIC':
            # If full frame data is needed, CDPs designed for a specific
            # subarray will be avoided.
            for sarray in MIRI_SUBARRAYS:
                avoid_strings.append(sarray)
            
        # The file name always contains the CDP type. A "_" is appended
        # to the match string to prevent cdptype "MASK" from matching the
        # coronographic MASK filters.
        if cdptype is not None and cdptype != 'ANY':
            match_strings.append(cdptype + "_")
                
        # An integration number is optional, but is always explicitly
        # included.  
        if integration is not None:
            match_strings.append(str(integration))
            
        # Find a list of files matching the given criteria except
        # the version number. match_cdp_substrings returns s sorted list.
        matched_files = self.match_cdp_substrings(mustcontain=match_strings,
                                                  mustnotcontain=avoid_strings)
        self.logger.debug("Matched files: " + str(matched_files))

        # If more than one file has been matched, a band has been specified
        # but a channel has not been provided, filter the strings to ensure
        # the band is preceeded by any numerical channel number or "_".
        # This is done to prevent a band string (such as LONG) being
        # accidentally matched against a detector name string (such as
        # MIRIFULONG). The test will not match MIRIFULONG but will match
        # 34LONG, 3LONG, _LONG, etc...
        if (len(matched_files) > 1) and \
           (band is not None and band != 'ANY') and \
           (channel is None or channel == 'ANY'):
            regexp_string = r"[A-Za-z0-9_\-]*[0-9_]" + band + "[A-Za-z0-9_\-]*"

            # Filter the list of files with this regular expression
            self.logger.debug("Version filtering with regexp: " + regexp_string)
            filtered_list = self._filter_regexp(matched_files, regexp_string,
                                                flags=re.IGNORECASE)
            matched_files = filtered_list

        # The list must now be sorted again to ensure that the CDP release
        # numbers are returned in the correct order (with the full release
        # versions without the letter suffix having a higher priority than the
        # prerelease versions ending in a letter).
        if len(matched_files) > 1:
            # This first set of releases will be valid from CDP-1 to CDP-9
            dig_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
            chr_list = ['A', 'B', 'C']
            release_order = []
            for dig in dig_list:
                for chr in chr_list:
                    release_order.append( '_' + dig + chr + '.' )
                release_order.append( '_' + '0' + dig + '.' )
            # This list extends the test up to CDP-10 to CDP-19, just in case.
            # If the releases go further than this, just add to the ten_list.
            ten_list = ['1']
            dig_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
            chr_list = ['A', 'B', 'C']
            for ten in ten_list:
                for dig in dig_list:
                    for chr in chr_list:
                        release_order.append( '_' + ten + dig + chr + '.' )
                    release_order.append( '_' + ten + dig + '.' )
                
            sorted_by_release = []
            list_copy = copy.copy(matched_files)
            for release_test in release_order:
                for filename in list_copy:
                    if release_test in filename:
                        # Add each file to the new list and remove it from the old
                        # list
                        sorted_by_release.append(filename)
                        matched_files.remove(filename)
            # To ensure that no files are accidentally missed, the new list
            # is appended to the remnant of the old list.
            matched_files.extend( sorted_by_release )

        # If a required version number is provided, it needs to be filtered
        # using a regular expression.
        if cdprelease is None and cdpversion is None and \
           cdpsubversion is None:
            # No version restrictions specified. Return the full list.
            return matched_files
        else:
            # The string starts with "MIRI_" and then contains any
            # number of alphanumeric characters and underscores.
            regexp_string = r"^MIRI_[A-Za-z0-9_\-]*"
            if cdprelease is not None:
                cdprelease = self._encode_version(cdprelease)
                regexp_string += cdprelease + r"\."
            else:
                # The CDP release number is a 1 or 2 digit number which can be
                # followed by zero or 1 uppercase alphabetic characters
                # (which signifies a special release, e.g. B for beta)
                # followed by a '.'
                regexp_string += r'[0-9]{1,2}[A-Z]?\.'
            if cdpversion is not None:
                cdpversion = self._encode_version(cdpversion)
                regexp_string += cdpversion + r"\."
            else:
                # A CDP version is a 2 digit number followed by a '.'
                regexp_string += r'[0-9]{2}\.'
            if cdpsubversion is not None:
                cdpsubversion = self._encode_version(cdpsubversion)
                regexp_string += cdpsubversion
            else:
                # A CDP version is a 2 digit number followed by a '.'
                regexp_string += r'[0-9]{2}'
            # The version code is the last part of the file name before
            # the .fits.
            regexp_string += "\.fits$"

            # Filter the list of files with this regular expression
            self.logger.debug("Version filtering with regexp: " + regexp_string)
            filtered_list = self._filter_regexp(matched_files, regexp_string,
                                                flags=re.IGNORECASE)
            return filtered_list

    def match_cdp_latest(self, cdptype, model='FM', detector=None,
                         readpatt=None, channel=None, band=None,
                         mirifilter=None, subarray='FULL', integration=None,
                         cdprelease=None, cdpversion=None, cdpsubversion=None):
        """
        
        Match the most recent version of a CDP contained in the FTP folder
        matching the given criteria.
        
        NOTE: If several files match the given criteria, only one file name
        will be returned - the one at the end of the sorted list. It is
        recommended that all the required parameters are specified except
        the version numbers.

        :Parameters:
    
        cdptype: str
            Data type of CDP to be obtained. Must be one of the defined
            MIRI CDP types.
        model: str (optional)
            The MIRI model required ('VM', 'FM' or 'JPL').
            Defaults to 'FM'.
        detector: str (optional)
            The MIRI detector required ('MIRIMAGE', 'MIRIFUSHORT',
            'MIRIFULONG' or 'ANY'). The old names ('IM', 'SW' or 'LW')
            may be used to find CDPs prior to the CDP-3 release.
            By default, matches any detector.
        readpatt: str (optional)
            The MIRI readout pattern required ('FAST', 'SLOW' or 'ANY').
            By default, matches any readout pattern.  
        channel: str (optional)
            The MIRI MRS channel required ('1', '2', '3', '4', '12',
            '34' or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches any channel.  
        band: str (optional)
            The MIRI MRS band required ('SHORT', 'MEDIUM', 'LONG',
            'SHORT-MEDIUM', 'SHORT-LONG', 'MEDIUM-SHORT', etc..., or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches any band.  
        mirifilter: str (optional)
            The MIRI filter required ('F560W', 'F770W', 'F1000W', 'F1130W',
            'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W', 'F2550WR',
            'F1065C', 'F1140C', 'F1550C', 'F2300C', 'P750L', 'FLENS', 'FND',
            'OPAQUE', 'ANY' or 'GENERIC').
            Valid only for imager data, when detector is 'IM'.
            By default, matches any filter.  
            Explicitly specifying 'GENERIC' will only match CDPs which do not
            specify any filter and therefore work with any filter (as opposed
            to 'ANY', which can match a CDP which is specific to one filter only).
        subarray: str (optional)
            The MIRI subarray required ('FULL', 'MASK1140', 'MASK1550',
            'MASK1065', 'MASKLYOT', 'BRIGHTSKY', 'SUB256', 'SUB128',
            'SUB64', 'SLITLESSPRISM').
            Defaults to 'FULL', which will match full frame data.
        integration: int (optional)
            When a CDP is split between separate files for each integration,
            this parameter can be used to specify a particular integration.
            Valid for 'DARK' or 'LINEARITY' data only.
            By default matches any integration.
        cdprelease: int or str (optional)
            A string or integer describing the CDP release number required.
            If not specified, the latest release will be matched.
        cdpversion: int or str (optional)
            A string or integer describing the CDP format version number
            required.
            If not specified, the latest version will be matched.
        cdpsubversion: int or str (optional)
            A string or integer describing the CDP subversion number required.
            If not specified, the latest subversion will be matched.
        
        :Returns:
        
        filename: list str
            The name of the latest version of the CDP matching the given
            criteria. If no candidates are matched returns an empty string.
         
        """
        candidates = self.match_cdp_filename(cdptype, model=model,
                                             detector=detector,
                                             readpatt=readpatt,
                                             channel=channel, band=band,
                                             mirifilter=mirifilter,
                                             subarray=subarray,
                                             integration=integration,
                                             cdprelease=cdprelease,
                                             cdpversion=cdpversion,
                                             cdpsubversion=cdpsubversion)
        # If only one candidate is returned, that is the one required.
        # Otherwise sort the candidates and return the last one.
        if candidates:
            ncandidates = len(candidates)
            if ncandidates <= 1:
                return candidates[0]
            else:
                # Give a warning in case the wrong candidate is chosen.
                strg = " There are %d matching CDP candidates:\n" % ncandidates
                for candidate in candidates:
                    strg += "  \'%s\'\n" % candidate
                strg += "of which only one (the last in this list) will be selected."
                # A warning is only necessary when different kinds of CDP are included in the sort.
                # If the CDP names only differ in the version number, there will be no more than
                # 2 characters different on average.
                if _diff_count_list(candidates) > 2:
                    self.logger.warn(strg)
                else:
                    self.logger.debug(strg)
                return candidates[-1]
        else:
            # No candidates are matched.
            return ''
        
    def __str__(self):
        """
        
        Return a string describing the object
        
        """
        if self.sftp is not None:
            strg = "  Files are obtained from the SFTP folder \'%s\'\n" % self.ftp_path
            strg += "  and cached to local folder \'%s\' " % self.cdp_dir
        else:
            strg = "  Relying on local folder \'%s\' " % self.cdp_dir
 
        abspath = os.path.abspath(self.cdp_dir)
        if abspath != self.cdp_dir:
            strg += "= \'%s\'.\n" % abspath
        else:
            strg += ".\n"
        strg += "  %s products are available:\n" % len(self.cdp_files_available)
        for cdpfile in self.cdp_files_available:
            strg += "    " + cdpfile + "\n"
        strg += "  %s documents are available:\n" % len(self.cdp_docs_available)
        for cdpfile in self.cdp_docs_available:
            strg += "    " + cdpfile + "\n"

        return strg


# Only one instance of this class is allowed to exist.
@six.add_metaclass(Singleton)
class MiriCDPInterface(object):
    """
    
    A class which manages the interface between the MIRI pipeline or
    simulator software and the MIRI Calibration Data Products.
    
    The MIRI CDPs are available in a repository which is accessible by
    SFTP from anywhere in the world. If a new CDP is required, it is copied
    from that repository.
    
    A cache of MIRI CDPs is maintained on the local machine, so that each
    CDP only needs to be copied from the repository once.
    
    NOTE: If this class is used from a location which cannot access the
    SFTP site (e.g. no internet connection or local access restrictions)
    the CDP files must first be copied manually into the local cache
    directory.
    
    The location of the local cache directory may be defined explicitly in
    the local_path parameter. If local_path is not defined, the local cache
    is found by querying environment variables. If the environment variable
    CDP_DIR is defined, this is assumed to point directly to the local cache.
    If CDP_DIR is not defined, the CDP cache is assumed to be contained in a
    'CDP' sub-directory of the location referenced by environment variable
    MIRI_DIR. If neither environment variable is defined, the current
    directory, '.', is used. The names of the environment variables are
    defined with the cdp_env_name and miri_env_name parameters.
    
    The class may be instructed to look for CDP files within a search path
    of more than one FTP folder. Each folder is managed internally by a
    MiriCDPFolder class.
    
    :Parameters:
    
    ftp_host: str (optional)
        The address of the machine hosting the MIRI SFTP repository.
        Defaults to the repository at Leuven.
        Setting this to 'LOCAL' will force MIRI CDPs to be obtained from
        the local cache.
    ftp_path: str (optional)
        The path to the folder (or folders) on the SFTP host where the
        MIRI CDPs are held. A search path consisting of more than one
        folder may be specified, separated by a ":" delimiter. These
        folders will be checked in the order given, and the first CDP
        matching the given criteria used.
        Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
        Defaults to the default CDP repository at Leuven.
    ftp_user: str (optional)
        A username with which to access the ftp site.
        Defaults to 'miri'.
    ftp_passwd: str (optional)
        A password with which to access the ftp site.
        Defaults to prompting for the password.
    timeout: float (optional)
        A new connection timeout in seconds.
        Defaults to 15 seconds.
    local_path: str (optional)
        A local file path giving the MIRI CDP folder.
        If not given, the local path will be obtained from the
        environment variable specified in the cdp_env_name parameter,
        or failing that from the miri_env_name parameter.
        If local_path is not defined, and none of the named environment
        variables are defined, the local path defaults to '.'.
        Examples: '.',  '/home/me/MIRI_DHAS/MPipeline/Cal".
    cdp_env_name: str (optional)
        The name of the environment variable defining the MIRI CDP path.
        Defaults to 'CDP_DIR'.
        If the named environment variable is not defined, the miri_env_name
        parameter is used to obtain the top-level MIRI environment name.
    miri_env_name: str (optional)
        The name of the environment variable containing the top level MIRI
        path. CDPs are assumed to be contained in a '/CDP' subdirectory.
        Defaults to 'MIRI_ENV'.  
        
    :Environment:
  
    cdp_env_name (CDP_DIR):
        Environment variable defining the MIRI CDP path.
        See the cdp_env_name parameter.
    miri_env_name (MIRI_ENV)
        Environment variable defining the top level MIRI path.
        See the miri_env_name parameter.
        
    """
    # Class variable defining the previous password
    PREVIOUS_PW = ''
    
    # Class variables defining the default SFTP parameters
    FTP_HOST_DEFAULT = 'www.miricle.org'        
    FTP_USER_DEFAULT = 'miri'
    FTP_PATH_DEFAULT  = 'CDP'
#     FTP_PATH_DEFAULT  = 'CDPSIM:CDP:CDPTMP'
    FTP_PATH_SEARCH = ':'
                                             
    def __init__(self, ftp_host=None, ftp_path=None, ftp_user='miri',
                 ftp_passwd='', timeout=15.0, local_path=None,
                 cdp_env_name='CDP_DIR', miri_env_name='MIRI_ENV',
                 logger=LOGGER):
        """
        
        Initialises the MiriCDPInterface class.
        
        Parameters: See class doc string.

        """
        # Get a Python logger object
        self.logger = logger.getChild(self.__class__.__name__)
        self.logger.setLevel(level=LOGGING_LEVEL)
        self.logger.debug("Class %s created for ftp_host=\'%s\', ftp_path=\'%s\'."  % \
                          (self.__class__.__name__, str(ftp_host), str(ftp_path)))

        self.sftp = None
        self.ftp_ok = False

        # Set up the SFTP/CDP connection environment
        self._setup(ftp_host, ftp_path, ftp_user, ftp_passwd, timeout,
                    local_path, cdp_env_name, miri_env_name)     

        # Open the connection
        self._open()        
        # Create a search list of FTP folders and initialise the list of
        # available CDPs within in each one. If the local cache is being
        # used, there is only one nominal FTP folder.
        self.cdp_folder_list = []
        if self.ftp_host != 'LOCAL' and self.ftp_ok:
            for ftpp in self.ftp_path.split(MiriCDPInterface.FTP_PATH_SEARCH):
                cdp_folder = MiriCDPFolder( self.sftp, ftpp, self.cdp_dir )
                cdp_folder.update_cdp_list()
                self.cdp_folder_list.append(cdp_folder)
        else:
            cdp_folder = MiriCDPFolder( self.sftp, '/', self.cdp_dir )
            cdp_folder.update_cdp_list()
            self.cdp_folder_list.append(cdp_folder)
#         # Close the connection
#         self._close()

    def __del__(self):
        """
        
        Tidies up the MiriCDPInterface class.
        
        """
        self.logger.debug("Tidying up class %s."  % self.__class__.__name__)
        for cdp_folder in self.cdp_folder_list:
            del cdp_folder
        del self.cdp_folder_list
        self._close()
 
    def refresh(self, ftp_host=None, ftp_path=None, ftp_user=None,
                ftp_passwd=None, timeout=15.0, local_path=None,
                cdp_env_name='CDP_DIR', miri_env_name='MIRI_ENV'):
        """
        
        Refreshes the SFTP and CDP environmental parameters if any have changed.
        
        Parameters: See class doc string.

        """
        if ftp_host is None:
            self.ftp_host = MiriCDPInterface.FTP_HOST_DEFAULT
        if ftp_user is None:
            self.ftp_user = MiriCDPInterface.FTP_USER_DEFAULT
        if ftp_path is None:
            self.ftp_path = MiriCDPInterface.FTP_PATH_DEFAULT

        # Only execute this method if any of the SFTP or local cache
        # parameters have changed.
        if (ftp_host != self.ftp_host) or \
           ((ftp_path is not None) and (ftp_path != self.ftp_path)) or \
           ((timeout is not None) and (timeout != self.timeout)) or \
           ((local_path is not None) and (local_path != self.local_path)) or \
           ((cdp_env_name is not None) and \
            (cdp_env_name != self.cdp_env_name)) or \
           ((miri_env_name is not None) and \
            (miri_env_name != self.miri_env_name)) or \
           (ftp_user != self.ftp_user) or \
           ((ftp_passwd is not None) and (ftp_passwd != self.ftp_passwd)):

            # Close any open connection and remove the old CDP folder list.
            self._close()
            for cdp_folder in self.cdp_folder_list:
                del cdp_folder
            self.cdp_folder_list = []
                    
            # Set up the SFTP/CDP connection environment
            self.logger.debug("Class %s refreshed with ftp_host=\'%s\', ftp_path=\'%s\'."  % \
                              (self.__class__.__name__, str(ftp_host), str(ftp_path)))
            self._setup(ftp_host, ftp_path, ftp_user, ftp_passwd, timeout,
                        local_path, cdp_env_name, miri_env_name)     

            # Open the connection
            self._open()
            # Create a search list of FTP folders and initialise the list of
            # available CDPs within in each one. If the local cache is being
            # used, there is only one nominal FTP folder.
            if self.ftp_host != 'LOCAL' and self.ftp_ok:
                for ftpp in self.ftp_path.split(MiriCDPInterface.FTP_PATH_SEARCH):
                    cdp_folder = MiriCDPFolder( self.sftp, ftpp, self.cdp_dir )
                    cdp_folder.update_cdp_list()
                    self.cdp_folder_list.append(cdp_folder)
            else:
                cdp_folder = MiriCDPFolder( self.sftp, '/', self.cdp_dir )
                cdp_folder.update_cdp_list()
                self.cdp_folder_list.append(cdp_folder)
#             # Close the connection
#             self._close()
            
    def _setup(self, ftp_host, ftp_path, ftp_user, ftp_passwd, timeout,
               local_path, cdp_env_name, miri_env_name):
        """
        
        Helper function which establishes the SFTP and CDP environment
        
        Parameters: See class doc string.
        
        """
        if ftp_user is not None:
            self.ftp_user = str(ftp_user)
        else:
            self.ftp_user = MiriCDPInterface.FTP_USER_DEFAULT
        if ftp_host is not None:
            self.ftp_host = str(ftp_host)
        else:
            self.ftp_host = MiriCDPInterface.FTP_HOST_DEFAULT
        # Prompt for a password if one has not been provided
        if self.ftp_host != 'LOCAL' and not ftp_passwd:
            if MiriCDPInterface.PREVIOUS_PW:
                self.logger.debug("Remembering previous ftp password.")
                ftp_passwd = MiriCDPInterface.PREVIOUS_PW
            else:
                print("")
                while not ftp_passwd:
                    # getpass does not work under Windows.
                    if sys.platform.startswith("win"):
                        ftp_passwd = raw_input( \
                                    'Please enter ftp password for %s:' % \
                                    ftp_user)
                    else:
                        ftp_passwd = getpass.getpass( \
                                    'Please enter ftp password for %s: ' % \
                                    ftp_user)
                # A password has been obtained from the user 
                self.user_pwd = True
        else:
            # A fixed password has been provided 
            self.user_pwd = False

        # Remember the password for next time.
        MiriCDPInterface.PREVIOUS_PW = ftp_passwd
        
        self.ftp_passwd = ftp_passwd
        if ftp_path is not None:
            self.ftp_path = str(ftp_path)
        else:
            self.ftp_path = MiriCDPInterface.FTP_PATH_DEFAULT  
        self.ftp_ok = False
        self.timeout = timeout
               
        if cdp_env_name:
            self.cdp_env_name = str(cdp_env_name)
        if miri_env_name:
            self.miri_env_name = str(miri_env_name)
        if local_path:
            self.logger.debug("Local path specified: " + str(local_path))
            if local_path and not local_path.isspace():
                # A local path has been specified explicitly.
                self.local_path = str(local_path)
            else:
                # An empty or blank local path translates to the current directory
                self.local_path = '.'
            self.cdp_dir = self.local_path
        else:
            # First try the MIRI CDP environment variable.
            self.local_path = os.getenv(self.cdp_env_name, '')
            if self.local_path:
                self.logger.debug("CDP path specified: " + str(self.local_path))
                self.cdp_dir = self.local_path
            else:
                # If the CDP environment variable could not be translated,
                # try the MIRI environment variable, falling back
                # on the current directory as a default.
                self.local_path = os.getenv(miri_env_name, '.')
                self.logger.debug("MIRI path specified: " + str(self.local_path))
                # The CDP directory is a '/CDP' subdirectory within that
                # top-level directory.   
                self.cdp_dir = os.path.join(self.local_path, 'CDP')

    def _open(self):
        """
        
        Helper function which opens the sftp connection.
        
        It assumes all the connection parameters are already contained in
        class attributes.
        
        The function is compatible with pysftp 0.2.8 (without host key
        checking) and pysftp 0.2.9 (with hostkey checking). A missing
        hostkey is converted from a major error into a warning.
                
        """
        self.logger.debug("Opening connection")
        if self.ftp_host == 'LOCAL':
            self.stfp = None
            return

        self.ftp_ok = False
        try:
            # First attempt to connect with hostkey checking
            self.sftp = pysftp.Connection(self.ftp_host,
                                          username=self.ftp_user,
                                          password=self.ftp_passwd)
            # If you get this far without an exception, the SFTP connection
            # has been successful.
            self.ftp_ok = True
        except SSHException as e:
            # Check for a "hostkey not found" exception
            emsg = str(e)
            if 'hostkey' in emsg:
                # Disable hostkey checking, if available (pysftp V0.2.9 onwards)
                if hasattr(pysftp, 'CnOpts'):
                    strg = 20 * "*" + "\n"
                    strg += "*** No hostkey found for %s." % self.ftp_host
                    strg += " Trying again without a hostkey.\n"
                    strg += "*** To improve security in future runs, please enter\n"
                    strg += "*** \t ssh %s \n" % self.ftp_host
                    strg += "*** on the command line and answer the prompt"
                    strg += " to create a hostkey.\n"
                    strg += 60 * "*"
                    self.logger.warn(strg)
                    cnopts = pysftp.CnOpts()
                    cnopts.hostkeys = None
                    # Connect without host checking
                    self.sftp = pysftp.Connection(self.ftp_host,
                                             username=self.ftp_user,
                                             password=self.ftp_passwd,
                                             cnopts=cnopts)
                    # If you get this far without an exception, the SFTP connection
                    # has been successful.
                    self.ftp_ok = True
                else:
                    strg = "SSHException: SFTP connection failed: %s" % emsg
                    strg += "  Falling back to LOCAL host."
                    self.logger.error(strg)
                    self.sftp = None
            else:
                strg = "SSHException: SFTP connection failed: %s" % emsg
                strg += "  Falling back to LOCAL host."
                self.logger.error(strg)
                self.sftp = None
        except Exception as e:
            strg = "%s: SFTP connection failed: %s" % \
                (e.__class__.__name__, str(e))
            strg += "  Falling back to LOCAL host."
            self.logger.error(strg)
            self.sftp = None
            
        # Change the timeout, if one has been provided.
        if self.timeout is not None:
            self.sftp.timeout = self.timeout

    def _close(self):
        """
        
        Helper function which closes the sftp connection.
        
        It assumes all the connection parameters are already contained in
        class attributes.
                
        """
        self.logger.debug("Closing connection")
        if self.sftp is not None:
            self.sftp.close()
            self.sftp = None
               
    def match_cdp_regexp(self, match_string, flags=re.IGNORECASE):
        """
        
        Return a list, concatenated from all the FTP folders in the search
        path, of all the CDP files matching a given regular expression.
        
        :Parameters:
        
        match_string: raw string
            A regular expression to be matched.
            For example, r'[A-Za-z_]*Bad[A-Za-z_]*.+fits' matches
            all bad pixel mask file names.
        flags: int (optional)
            Flags used to modify the regular expression search.
            For example, the default flags=re.IGNORECASE will treat 
            lowercase and uppercase matches equally.
            Specific flags=0 to match the exact case.
            
        :Returns:
        
        (filenames: list of str, folders: list of str)
            A sorted list of matching filenames and corresponding folders
        
        """
        self.logger.debug("Matching regexp: " + match_string)
        # Make a concatenated list of all the files found in all
        # folders of the search path.
        matched_files = []
        matched_folders = []
        for cdp_folder in self.cdp_folder_list:
            files_subset = cdp_folder.match_cdp_regexp(match_string,
                                                       flags=flags)
            if files_subset:
                matched_files += files_subset
                matched_folders += [cdp_folder.ftp_path] * len(files_subset)
        return (matched_files, matched_folders)

    def match_cdp_substrings(self, mustcontain=[], mustnotcontain=[]):
        """
        
        Return a list, concatenated from all the FTP folders in the search
        path, of all the CDP files whose names contain all the compulsory
        substrings and do not contain any forbidden substrings.
        
        :Parameters:
        
        mustcontain: list of str
            A list of substrings which must be matched.
        mustnotcontain: list of str
            A list of substrings which must NOT be matched.
            
        :Returns:
        
        (filenames: list of str, folders: list of str)
            A sorted list of matching filenames and corresponding folders.
            Returns an empty filename list if no files are matched.
        
        """
        self.logger.debug("Must contain: " + str(mustcontain))
        self.logger.debug("Must not contain: " + str(mustnotcontain))
        # Make a concatenated list of all the files found in all
        # folders of the search path.
        matched_files = []
        matched_folders = []
        for cdp_folder in self.cdp_folder_list:
            files_subset = cdp_folder.match_cdp_substrings(
                                        mustcontain=mustcontain,
                                        mustnotcontain=mustnotcontain)
            if files_subset:
                matched_files += files_subset
                matched_folders += [cdp_folder.ftp_path] * len(files_subset)
        return (matched_files, matched_folders)

    def match_cdp_filename(self, cdptype, model='FM', detector=None,
                           readpatt=None, channel=None, band=None,
                           mirifilter=None, subarray='FULL', integration=None,
                           cdprelease=None, cdpversion=None, cdpsubversion=None):
        """
        
        Find a list, concatenated from all the FTP folders in the search
        path, of available CDPs matching the required criteria.
        The CDP file naming convention is assumed to be:
        
        MIRI_<model>_<detector>_<detsetng>_<readpatt>_<channelband_or_filter>_<subarray>_<reftype>_<version>.fits
        
        :Parameters:
    
        cdptype: str
            Data type of CDP to be obtained. Must be one of the defined
            MIRI CDP types. For example, 'MASK', 'DARK', 'LASTFRAME',
            'DISTORTION', 'DROOP', 'PIXELFLAT', 'FRINGE', 'PHOTOM', 'IPC',
            'COLCORR', 'FLUX', 'JUMP', 'LATENT', 'LINEARITY', 'SATURATION',
            'PSF', 'STRAY', 'TRACORR', 'WAVCORR', 'TELEM', etc....
        model: str (optional)
            The MIRI model required ('VM', 'FM' or 'JPL').
            Defaults to 'FM'.
        detector: str (optional)
            The MIRI detector required ('MIRIMAGE', 'MIRIFUSHORT',
            'MIRIFULONG' or 'ANY'). The old names ('IM', 'SW' or 'LW')
            may be used to find CDPs prior to the CDP-3 release.
            By default, matches all detectors.
        readpatt: str (optional)
            The MIRI readout pattern required ('FAST', 'SLOW' or 'ANY').
            By default, matches all readout patterns.  
        channel: str (optional)
            The MIRI MRS channel required ('1', '2', '3', '4', '12',
            '34' or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches all channels.  
        band: str (optional)
            The MIRI MRS band required ('SHORT', 'MEDIUM', 'LONG',
            'SHORT-MEDIUM', 'SHORT-LONG', 'MEDIUM-SHORT', etc..., or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches all bands.  
        mirifilter: str (optional)
            The MIRI filter required ('F560W', 'F770W', 'F1000W', 'F1130W',
            'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W', 'F2550WR',
            'F1065C', 'F1140C', 'F1550C', 'F2300C', 'P750L', 'FLENS', 'FND',
            'OPAQUE', 'ANY' or 'GENERIC').
            Valid only for imager data, when detector is 'IM'.
            By default, matches all filters.  
            Explicitly specifying 'GENERIC' will only match CDPs which do not
            specify any filter and therefore work with any filter (as opposed
            to 'ANY', which can match a CDP which is specific to one filter only).
        subarray: str (optional)
            The MIRI subarray required ('FULL', 'MASK1140', 'MASK1550',
            'MASK1065', 'MASKLYOT', 'BRIGHTSKY', 'SUB256', 'SUB128',
            'SUB64', 'SLITLESSPRISM').
            Defaults to 'FULL', which will match full frame data.
        integration: int (optional)
            When a CDP is split between separate files for each integration,
            this parameter can be used to specify a particular integration.
            Valid for 'Dark' data only.
            By default matches all integrations.
        cdprelease: int or str (optional)
            A string or integer describing the CDP release number required.
            If not specified, the latest release will be matched.
        cdpversion: int or str (optional)
            A string or integer describing the CDP format version number
            required.
            If not specified, the latest version will be matched.
        cdpsubversion: int or str (optional)
            A string or integer describing the CDP subversion number required.
            If not specified, the latest subversion will be matched.
        
        :Returns:
        
        (filenames: list of str, folders: list of str)
            A sorted list of names of files, and their corresponding folders,
            matching the given criteria.
            Returns an empty filename list if no files are matched.
         
        """
        # Make a concatenated list of all the files found in all
        # folders of the search path.
        matched_files = []
        matched_folders = []
        for cdp_folder in self.cdp_folder_list:
            files_subset = cdp_folder.match_cdp_filename(cdptype, model=model,
                                        detector=detector, readpatt=readpatt,
                                        channel=channel, band=band,
                                        mirifilter=mirifilter,
                                        subarray=subarray,
                                        integration=integration,
                                        cdprelease=cdprelease,
                                        cdpversion=cdpversion,
                                        cdpsubversion=cdpsubversion)
            if files_subset:
                matched_files += files_subset
                matched_folders += [cdp_folder.ftp_path] * len(files_subset)
        return (matched_files, matched_folders)

    def match_cdp_latest(self, cdptype, model='FM', detector=None,
                         readpatt=None, channel=None, band=None,
                         mirifilter=None, subarray='FULL', integration=None,
                         cdprelease=None, cdpversion=None, cdpsubversion=None):
        """
        
        Match the most recent version of a CDP matching the given criteria,
        searching all the FTP folders in the given search path.
        
        NOTE: If several files match the given criteria, only one file name
        will be returned - the one at the end of the sorted list. It is
        recommended that all the required parameters are specified except
        the version numbers.

        :Parameters:
    
        cdptype: str
            Data type of CDP to be obtained. Must be one of the defined
            MIRI CDP types.
        model: str (optional)
            The MIRI model required ('VM', 'FM' or 'JPL').
            Defaults to 'FM'.
        detector: str (optional)
            The MIRI detector required ('MIRIMAGE', 'MIRIFUSHORT',
            'MIRIFULONG' or 'ANY'). The old names ('IM', 'SW' or 'LW')
            may be used to find CDPs prior to the CDP-3 release.
            By default, matches any detector.
        readpatt: str (optional)
            The MIRI readout pattern required ('FAST', 'SLOW' or 'ANY').
            By default, matches any readout pattern.  
        channel: str (optional)
            The MIRI MRS channel required ('1', '2', '3', '4', '12',
            '34' or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches any channel.  
        band: str (optional)
            The MIRI MRS band required ('SHORT', 'MEDIUM', 'LONG',
            'SHORT-MEDIUM', 'SHORT-LONG', 'MEDIUM-SHORT', etc..., or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches any band.  
        mirifilter: str (optional)
            The MIRI filter required ('F560W', 'F770W', 'F1000W', 'F1130W',
            'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W', 'F2550WR',
            'F1065C', 'F1140C', 'F1550C', 'F2300C', 'P750L', 'FLENS', 'FND',
            'OPAQUE', 'ANY' or 'GENERIC').
            Valid only for imager data, when detector is 'IM'.
            By default, matches any filter.  
            Explicitly specifying 'GENERIC' will only match CDPs which do not
            specify any filter and therefore work with any filter (as opposed
            to 'ANY', which can match a CDP which is specific to one filter only).
        subarray: str (optional)
            The MIRI subarray required ('FULL', 'MASK1140', 'MASK1550',
            'MASK1065', 'MASKLYOT', 'BRIGHTSKY', 'SUB256', 'SUB128',
            'SUB64', 'SLITLESSPRISM').
            Defaults to 'FULL', which will match full frame data.
        integration: int (optional)
            When a CDP is split between separate files for each integration,
            this parameter can be used to specify a particular integration.
            Valid for 'DARK' or 'LINEARITY' data only.
            By default matches any integration.
        cdprelease: int or str (optional)
            A string or integer describing the CDP release number required.
            If not specified, the latest release will be matched.
        cdpversion: int or str (optional)
            A string or integer describing the CDP format version number
            required.
            If not specified, the latest version will be matched.
        cdpsubversion: int or str (optional)
            A string or integer describing the CDP subversion number required.
            If not specified, the latest subversion will be matched.
        
        :Returns:
        
        (filename: str, ftp_path: str)
            The name of the latest version of the CDP matching the given
            criteria. If no candidates are matched returns an empty string.
         
        """
        for cdp_folder in self.cdp_folder_list:
            candidate = cdp_folder.match_cdp_latest(cdptype,
                            model=model, detector=detector,
                            readpatt=readpatt, channel=channel, band=band,
                            mirifilter=mirifilter, subarray=subarray,
                            integration=integration,
                            cdprelease=cdprelease, cdpversion=cdpversion,
                            cdpsubversion=cdpsubversion)
            if candidate:
                return (candidate, cdp_folder.ftp_path)
        # The search has failed.
        return ('', '')


    def update_cache(self, cdp_files='', ftp_path='/', force_update=False):
        """
        
        Bring the cached version of the specified file up to date by copying
        it over from the SFTP site.
        
        NOTE: How do you know when a file on the SFTP site is newer than
        the one stored in the local cache?
        
        :Parameters:
        
        cdp_files: str or list of str
            Either a single file name or a list of names of files to
            be updated. If an empty string is provided ALL available
            CDP files will be copied to the cache.
        ftp_path: str (optional)
            A file path on the FTP repository.
            Defaults to the top-level folder, '/'.
        force_update: bool (optional)
            Files will only be pulled across from the SFTP repository
            if they don't exist in the cache. Setting force_update to
            True will force the files to be updated and overwritten
            regardless.
        
        :Returns:
        
        local_filename: str
            The name of the last file to be matched or updated. Empty
            string if no files were matched.
         
        """
        # If a filename or list of files is not specified, copy all known CDPs.
        if cdp_files is None or not cdp_files:
            cdp_files = []
            for cdp_folder in self.cdp_folder_list:
                cdp_files += cdp_folder.cdp_files_available
        elif not isinstance(cdp_files, (tuple,list)):
            # If a list has not been provided, make one
            cdp_files = [str(cdp_files)]

        # Create the local path if it doesn't already exist.
        if not os.path.isdir(self.cdp_dir):
            self.logger.debug("Making new directory, %s", self.cdp_dir)
            os.makedirs(self.cdp_dir)
                    
        local_filename = ''
#         open_attempted = False
        for cdp_file in cdp_files:
        
            if cdp_file:
                local_filename = os.path.join(self.cdp_dir, cdp_file)
                # Check if the local file already exists.
                # TODO: Check if the local file is older than the one in the repository.
                if not os.path.isfile(local_filename) or force_update:
#                     # Open the connection the first time an attempt is made to
#                     # access a file which isn't in the local cache.
#                     if not open_attempted:
#                         self._open()
#                         open_attempted = False
                    
                    # The local path needs translating into something SFTP understands.
                    public_url = self.ftp_host + '/.../' + cdp_file
                    strg = " Retrieving CDP file\n  from sftp \'%s\'\n" % public_url
                    strg += "  to cache \'%s\' ..." % local_filename
                    self.logger.info(strg)
                    if self.sftp is not None:
                        new_local_filename = local_filename.replace(os.path.sep, '/')
                        try:
                            self.sftp.chdir(ftp_path)
                        except IOError, e:
                            strg = "IOError: Failed to change directory to FTP folder \'%s\'\n" % ftp_path
                            strg += "  %s" % str(e)
                            raise IOError(strg)
                        self.sftp.get(cdp_file, localpath=new_local_filename)
                        self.sftp.chdir('/')
                        # TODO: Check the integrity of the copied file using a checksum.
                    else:
                        strg = " Cached CDP file \'%s\' has been removed!" % local_filename
                        strg += " It cannot be retrieved from LOCAL."
                        self.logger.error(strg)
                else:
                    if os.stat(local_filename).st_size > 0:
                        strg = " Cached CDP file \'%s\' already exists." % local_filename
                        self.logger.debug(strg)
                    else:
                        strg = " Cached CDP file \'%s\' is empty!" % local_filename
                        strg += " Delete this file and try again."
                        self.logger.error(strg)
            else:
                strg = " Ignoring empty file name \'%s\'" % str(cdp_file)
                self.logger.warn(strg)
#         # Ensure the connection is closed.
#         self._close()
        return local_filename

    def get_cdp_file(self, cdptype, model='FM', detector=None, readpatt=None,
                     channel=None, band=None, mirifilter=None, subarray='FULL',
                     integration=None, cdprelease=None, cdpversion=None,
                     cdpsubversion=None, force_update=False):
        """
        
        Get a specified CDP, copy it to the local cache if necessary
        and return its file path so it may be opened.

        NOTE: If several files match the given criteria, only one file
        will be returned.
        
        :Parameters:
    
        cdptype: str
            Data type of CDP to be obtained. Must be one of the defined
            MIRI CDP types. For example, 'MASK', 'DARK', 'LASTFRAME',
            'DISTORTION', 'DROOP', 'PIXELFLAT', 'FRINGE', 'PHOTOM', 'IPC',
            'COLCORR', 'FLUX', 'JUMP', 'LATENT', 'LINEARITY', 'SATURATION',
            'PSF', 'STRAY', 'TRACORR', 'WAVCORR', 'TELEM', etc....
        model: str (optional)
            The MIRI model required ('VM', 'FM' or 'JPL').
            Defaults to 'FM'.
        detector: str (optional)
            The MIRI detector required ('MIRIMAGE', 'MIRIFUSHORT',
            'MIRIFULONG' or 'ANY'). The old names ('IM', 'SQ' or 'LW')
            may be used to find CDPs prior to the CDP-3 release.
            By default, matches all detectors.
        readpatt: str (optional)
            The MIRI readout pattern required ('FAST', 'SLOW' or 'ANY').
            By default, matches all readout patterns.  
        channel: str (optional)
            The MIRI MRS channel required ('1', '2', '3', '4' or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches all channels.  
        band: str (optional)
            The MIRI MRS band required ('A', 'B', 'C', or 'ANY').
            Valid only for MRS data, when detector is 'MIRIFUSHORT' or
            'MIRIFULONG'.
            By default, matches all bands.  
        mirifilter: str (optional)
            The MIRI filter required ('F560W', 'F770W', 'F1000W', 'F1130W',
            'F1280W', 'F1500W', 'F1800W', 'F2100W', 'F2550W', 'F2550WR',
            'F1065C', 'F1140C', 'F1550C', 'F2300C', 'P750L', 'FLENS', 'FND',
            'OPAQUE', 'ANY' or 'GENERIC').
            Valid only for imager data, when detector is 'IM'.
            By default, matches all filters.  
            Explicitly specifying 'GENERIC' will only match CDPs which do not
            specify any filter and therefore work with any filter (as opposed
            to 'ANY', which can match a CDP which is specific to one filter only).
        subarray: str (optional)
            The MIRI subarray required ('FULL', 'MASK1140', 'MASK1550',
            'MASK1065', 'MASKLYOT', 'BRIGHTSKY', 'SUB256', 'SUB128',
            'SUB64', 'SLITLESSPRISM').
            Defaults to 'FULL', which will match full frame data..
        integration: int (optional)
            When a CDP is split between separate files for each integration,
            this parameter can be used to specify a particular integration.
            Valid for 'Dark' data only.
            By default matches all integrations.
        cdprelease: int or str (optional)
            A string or integer describing the CDP release number required.
            If not specified, the latest release will be matched.
        cdpversion: int or str (optional)
            A string or integer describing the CDP format version number
            required.
            If not specified, the latest version will be matched.
        cdpsubversion: int or str (optional)
            A string or integer describing the CDP subversion number required.
            If not specified, the latest subversion will be matched.
        force_update: bool (optional)
            Files will only be pulled across from the SFTP repository
            if they don't exist in the cache. Setting force_update to
            True will force the files to be updated and overwritten
            regardless.
        
        :Returns:
        
        local_filename: str
            The name of the local file created or matched. Empty
            string if no files were matched.
        
        """
        (matching_file, ftp_path) = self.match_cdp_latest(cdptype, model=model,
                                              detector=detector,
                                              readpatt=readpatt,
                                              channel=channel, band=band,
                                              mirifilter=mirifilter,
                                              subarray=subarray,
                                              integration=integration,
                                              cdprelease=cdprelease,
                                              cdpversion=cdpversion,
                                              cdpsubversion=cdpsubversion)
        if not matching_file:
            strg = "No CDP file matching the specified criteria."
            raise ValueError(strg)
        
        local_filename = self.update_cache(matching_file, ftp_path,
                                           force_update=force_update)
        return local_filename
    
    def get_cdp_doc(self, prefix, issue='', force_update=False):
        """
        
        Get one or more CDP documents matching the specified prefix,
        copying each file to the local cache if necessary and return
        a list of the files copied.
        
        Typically, only one document should match a given unique prefix,
        but there might be several issues of the same document or a set
        of documents on the same subject. The "issue" parameter can be used
        to narrow the search to a particular issue. If more than one
        document is returned, the reader should check each document to
        determine which is the most relevant.
        
        :Parameters:
    
        prefix: str
            Prefix of the document file name. For example a prefix of
            "MIRI-TN-00001-ETH" will match any file starting with this
            string.
        issue: str (optional)
            An optional issue string which must also be matched.
            The default is no string, which will match all issues.
        force_update: bool (optional)
            Files will only be pulled across from the SFTP repository
            if they don't exist in the cache. Setting force_update to
            True will force the files to be updated and overwritten
            regardless.
        
        :Returns:
        
        docs_copied: list of str
            The name(s) of the documents matched and copied to the cache.
            Empty list if no files were matched.
        
        """
        matched_files = []
        for cdp_folder in self.cdp_folder_list:
            for cdp_doc in cdp_folder.cdp_docs_available:
#                 self.logger.debug("Testing ftp_file=%s", ftp_file)
                if re.match("^" + prefix, cdp_doc) is not None:
                    if (len(issue) < 1) or (issue in cdp_doc):
                        matched_files.append(cdp_doc)                
            matched_files.sort()
            
            files_copied = []
            open_attempted = False
            for cdp_doc in matched_files:
            
                if cdp_doc:
                    local_filename = os.path.join(self.cdp_dir, cdp_doc)
                    files_copied.append(local_filename)
                    # Check if the local file already exists.
                    # TODO: Check if the local file is older than the one in the repository.
                    if not os.path.isfile(local_filename) or force_update:
                        # Open the connection the first time an attempt is made to
                        # access a file which isn't in the local cache.
                        if not open_attempted:
                            self._open()
                            open_attempted = False
                        
                        # The local path needs translating into something SFTP understands.
                        public_url = self.ftp_host + '/.../' + cdp_file
                        strg = " Retrieving CDP document\n  from sftp \'%s\'\n" % public_url
                        strg += "  to cache \'%s\' ..." % local_filename
                        self.logger.info(strg)
                        new_local_filename = local_filename.replace(os.path.sep, '/')
                        try:
                            self.sftp.chdir(cdp_folder.ftp_path)
                        except IOError, e:
                            strg = "IOError: Failed to change directory to FTP folder \'%s\'\n" % cdp_folder.ftp_path
                            strg += "  %s" % str(e)
                            raise IOError(strg)
                        self.sftp.chdir('/')
                        self.sftp.get(cdp_doc, localpath=new_local_filename)
                        # TODO: Check the integrity of the copied file using a checksum.
                    else:
                        strg = " Cached CDP document \'%s\' already exists." % local_filename
                        self.logger.debug(strg)
                else:
                    strg = " Ignoring empty document name \'%s\'" % str(cdp_doc)
                    self.logger.warn(strg)
        # Ensure the connection is closed.
        self._close()
        return files_copied
        
    def __str__(self):
        """
        
        Return a string describing the object
        
        """
        strg =  "MIRI Calibration Data Product Interface\n"
        strg += "=======================================\n"
        if self.ftp_ok and (self.sftp is not None):
            strg += "Files are obtained from the SFTP site \'%s\'\n" % self.ftp_host
            strg += "and cached to local folder \'%s\' " % self.cdp_dir
        elif self.ftp_host == 'LOCAL':
            strg += "LOCAL host. Relying on local folder \'%s\' " % self.cdp_dir
        else:
            strg += "Unable to connect to SFTP site \'%s\' !!\n" % self.ftp_host
            strg += "so relying on local folder \'%s\' " % self.cdp_dir

        abspath = os.path.abspath(self.cdp_dir)
        if abspath != self.cdp_dir:
            strg += "= \'%s\'.\n" % abspath
        else:
            strg += ".\n"

        nfolders = len(self.cdp_folder_list)
        if nfolders == 1:
            strg += "There is just 1 FTP folder in the search path.\n"
        else:
            strg += "There are %d FTP folders in the search path.\n" % nfolders
        seqnum = 1
        for cdp_folder in self.cdp_folder_list:
            strg += "SEARCH PATH (%d)\n" % seqnum
            strg += "----------------\n"
            strg += str(cdp_folder)
            seqnum += 1
        strg += "\n"

        return strg


#
# A minimal test and some examples of how to use the above utilities
# are run when this file is executed as a main program.
#
if __name__ == '__main__':
    import time    
    SHORT_DELAY = 0.2  # Short delay to let log output complete before printing        
    MEDIUM_DELAY = 1.0 # Medium delay to let large print output catch up before logging
    LONG_DELAY = 10.0  # Long delay to prevent the SFTP server from being flooded with requests       
    
    print("Testing the cdp utilities module.")
    
    # Set these variables to enable or disable some of the tests.
    TEST_MATCH = False   # Test the FTP query and string pattern matching
    TEST_GET = False     # Test the get_cdp function
    TEST_SIM = True      # Test the get_cdp function as used by MIRI simulators
    TEST_ERRORS = False  # Test the functions respond sensibly to errors
    # Set this variable False to turn off the getting of the files
    # and confine the testing just to pattern matching.
    GET_FILES = True
    
    # Set this variable False to suppress the plotting
    PLOTTING = False
    # Set these variable False to suppress detailed output
    VERBOSE_CDPIFACE = True  # Display contents of CDP interface objects
    VERBOSE_MODELS = False   # Display contents of CDP data objects

    if TEST_MATCH:
        print("\nTesting the MiriCDPInterface class and file matching functions")
        # Create a MIRI CDP interface object and print a summary.
        miricdp = MiriCDPInterface()
        if VERBOSE_CDPIFACE:
            print( miricdp )
            # Wait for the print output catch up.
            time.sleep(MEDIUM_DELAY)
           
        # Get a list of old bad pixel masks by matching a regular expression
        print( "\nOld bad pixel masks from regexp:" )
        (filenames, folders) = \
            miricdp.match_cdp_regexp(r'MIRI_[A-Za-z_]*Bad[A-Za-z_]*.+fits')
        print( str(filenames) )
        time.sleep(MEDIUM_DELAY)
    
        # Check availability of flat-fields.
        print( "\nFlat-fields available for various combinations" )
        for detector in ('MIRIMAGE', 'MIRIFUSHORT', 'MIRIFULONG'):
            for readpatt in ('FAST', 'SLOW', 'ANY'):
                for subarray in ('FULL', 'BRIGHTSKY'):
                    print(detector, "with READPATT=", readpatt,
                          "and subarray=", subarray)
                    (newestflat, ftp_path) = miricdp.match_cdp_latest('FLAT',
                                                detector=detector,
                                                readpatt=readpatt,
                                                subarray=subarray)
                    if newestflat:
                        print( newestflat, ftp_path )
                    else:
                        print( "<not found>")
        time.sleep(SHORT_DELAY)
       
        # Get a list of bad pixel masks by matching specified criteria.
        print( "\nFile names of new bad pixel masks for the imager:" )
        (filenames, folders) = \
            miricdp.match_cdp_filename('MASK', detector='MIRIMAGE')
        print( str(filenames) )
        time.sleep(SHORT_DELAY)
    
        # Get the most recent bad pixel mask matching the specified criteria.
        print( "\nMost recent bad pixel mask for the imager" )
        (newestbad, ftp_path) = miricdp.match_cdp_latest('MASK',
                                                         detector='MIRIMAGE')
        print( newestbad, ftp_path )
        time.sleep(SHORT_DELAY)
        if GET_FILES and newestbad:
            print( "Updating cache with this bad pixel mask..." )
            last_file = miricdp.update_cache(newestbad, ftp_path=ftp_path)
            time.sleep(LONG_DELAY)
    
        print( "\nMost recent release 3 bad pixel mask for the imager" )
        (newestbad, ftp_path) = miricdp.match_cdp_latest('MASK',
                                            detector='MIRIMAGE', cdprelease=3)
        print( newestbad, ftp_path )
        time.sleep(SHORT_DELAY)
        if GET_FILES and newestbad:
            print( "Updating cache with this bad pixel mask..." )
            last_file = miricdp.update_cache(newestbad, ftp_path=ftp_path)
            time.sleep(LONG_DELAY)
    
        print( "\nMost recent dark for the imager." )
        (newestdark, ftp_path) = miricdp.match_cdp_latest('DARK', detector='MIRIMAGE')
        print( newestdark, ftp_path )
        time.sleep(SHORT_DELAY)
     
        print( "\nMost recent dark for the imager for the FAST readout pattern with any subarray mode." )
        (newestdark, ftp_path) = miricdp.match_cdp_latest('DARK', detector='MIRIMAGE',
                                            readpatt='FAST', subarray='ANY')
        print( newestdark, ftp_path )
        time.sleep(SHORT_DELAY)
     
        print( "\nMost recent FULL frame dark for the imager for the FAST readout pattern." )
        (newestdark, ftp_path) = miricdp.match_cdp_latest('DARK', detector='MIRIMAGE',
                                            readpatt='FAST', subarray='FULL')
        print( newestdark, ftp_path )
        time.sleep(SHORT_DELAY)
     
        print( "\nMost recent dark for the imager for the MASK1150 subarray." )
        (newestdark, ftp_path) = miricdp.match_cdp_latest('DARK', detector='MIRIMAGE',
                                            subarray='MASK1550')
        print( newestdark, ftp_path )
        time.sleep(SHORT_DELAY)
     
        print( "\nMost recent pixel flat-field for the imager with F770W filter" )
        (newestflat, ftp_path) = miricdp.match_cdp_latest('FLAT', detector='MIRIMAGE',
                                            mirifilter='F770W')
        print( newestflat, ftp_path )
        time.sleep(SHORT_DELAY)
     
        print( "\nMost recent FULL frame pixel flat-field for the imager with F770W filter" )
        (newestflat, ftp_path) = miricdp.match_cdp_latest('FLAT', detector='MIRIMAGE',
                                            mirifilter='F770W', subarray='FULL' )
        print( newestflat, ftp_path )
        time.sleep(SHORT_DELAY)
    
        # Get a list of all the CDPs for the JPL model.
        print( "\nFile names for all the JPL model CDPs" )
        (jpl_list, folder_list) = miricdp.match_cdp_filename('ANY', model='JPL')
        print( str(jpl_list) )
        time.sleep(SHORT_DELAY)
    
        print( "\nFile names of any CDPs for the SUB64 subarray" )
        (sub64_list, folder_list) = miricdp.match_cdp_filename('ANY', detector='ANY',
                                            subarray='SUB64')
        print( str(sub64_list) )
        time.sleep(SHORT_DELAY)
         
        # Refresh the CDP object to use the local cache only
        miricdp.refresh(ftp_host='LOCAL')
        if VERBOSE_CDPIFACE:
            print( miricdp )
            # Wait for the print output catch up.
            time.sleep(MEDIUM_DELAY)
        del miricdp
    
    if TEST_GET and GET_FILES:
        print("\nTesting general use of get_cdp")
        print( "Getting the latest release 5 bad pixel mask for the imager from the CDP repository" )
        datamodel = get_cdp('MASK', detector='MIRIMAGE', cdprelease=5)
        if datamodel is not None:
            if VERBOSE_MODELS:
                print( datamodel )
            else:
                print( datamodel.__class__.__name__, " obtained successfully." )                
            if PLOTTING:
                datamodel.plot("Latest release 5 bad pixel mask")
            del datamodel
        time.sleep(LONG_DELAY)

        print( "Getting the latest pixel flats for the imager from the CDP repository" )
        datamodel = get_cdp('PIXELFLAT', detector='MIRIMAGE')
        if datamodel is not None:
            if VERBOSE_MODELS:
                print( datamodel )
            else:
                print( datamodel.__class__.__name__, " obtained successfully." )                
            if PLOTTING:
                datamodel.plot("Latest pixel flat")
            del datamodel
        time.sleep(LONG_DELAY)

        print( "Getting the gain model for the imager from the CDP repository" )
        datamodel = get_cdp('GAIN', detector='MIRIMAGE')
        if datamodel is not None:
            if VERBOSE_MODELS:
                print( datamodel )
            else:
                print( datamodel.__class__.__name__, " obtained successfully." )                
            if PLOTTING:
                datamodel.plot("Latest gain map")
            del datamodel
        time.sleep(LONG_DELAY)
  
        print( "Getting the readnoise model for the imager from the CDP repository" )
        datamodel = get_cdp('READNOISE', detector='MIRIMAGE')
        if datamodel is not None:
            if VERBOSE_MODELS:
                print( datamodel )
            else:
                print( datamodel.__class__.__name__, " obtained successfully." )                
            if PLOTTING:
                datamodel.plot("Latest readnoise map")
            del datamodel
        time.sleep(LONG_DELAY)
         
        print( "Getting the release 2 SRF for the LRS from the CDP repository" )
        datamodel = get_cdp('SRF', detector='IM', mirifilter='P750L')
        if datamodel is not None:
            if VERBOSE_MODELS:
                print( datamodel )
            else:
                print( datamodel.__class__.__name__, " obtained successfully." )                
            if PLOTTING:
                datamodel.plot("Latest LRS SRF")
            del datamodel
        time.sleep(LONG_DELAY)

    if TEST_SIM and GET_FILES:
        print("\nTesting use of get_cdp by the simulators")
        cdps_found = []
        cdps_not_found = []

        #ftp_path = 'CDPSIM'
        FTP_PATH = '/CDPSIM/1.0/:/CDPSIM/'
        FTP_USER = 'cdpuser'
        FTP_PASS = 'R7ZWEXEEsAH7'
        miricdp = MiriCDPInterface(ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                   ftp_passwd=FTP_PASS)
        if VERBOSE_CDPIFACE:
            print( miricdp )
            # Wait for the print output catch up.
            time.sleep(MEDIUM_DELAY)
        del miricdp

        print("Bad pixel masks")
        for detector in ('MIRIMAGE', 'MIRIFUSHORT', 'MIRIFULONG'):
            strg = "MASK for detector=%s" % detector
            print("\n" + strg)
            maskmodel = get_cdp('MASK', detector=detector, ftp_path=FTP_PATH,
                                ftp_user=FTP_USER, ftp_passwd=FTP_PASS)
            if maskmodel is not None:
                cdps_found.append(strg)
                if VERBOSE_MODELS:
                    print( maskmodel )
                if PLOTTING:
                    maskmodel.plot(strg)
            else:
                cdps_not_found.append(strg)
                print("*** CDP NOT FOUND ***")
            del maskmodel
            time.sleep(LONG_DELAY)
 
        print("Gain models")
        for detector in ('MIRIMAGE', 'MIRIFUSHORT', 'MIRIFULONG'):
            strg = "GAIN for detector=%s" % detector
            print("\n" + strg)
            gainmodel = get_cdp('GAIN', detector=detector, ftp_path=FTP_PATH,
                                ftp_user=FTP_USER, ftp_passwd=FTP_PASS)
            if gainmodel is not None:
                cdps_found.append(strg)
                if VERBOSE_MODELS:
                    print( gainmodel )
                if PLOTTING:
                    gainmodel.plot(strg)
            else:
                cdps_not_found.append(strg)
                print("*** CDP NOT FOUND ***")
            del gainmodel
            time.sleep(LONG_DELAY)

        def find_dark( detector, readpatt, subarray, averaged=False):
            """
            
            Helper function to find a DARK CDP in the same way as is
            done within SCASim.
            
            """
            from miri.datamodels.cdp import MiriDarkReferenceModel
            miricdp = MiriCDPInterface(ftp_path=FTP_PATH,
                                       ftp_user=FTP_USER,
                                       ftp_passwd=FTP_PASS)
            miricdp.refresh(ftp_path=FTP_PATH,
                            ftp_user=FTP_USER, 
                            ftp_passwd=FTP_PASS)
            must_contain = []
            must_not_contain = []
            if detector:
                must_contain.append(detector)
            if readpatt:
                must_contain.append(readpatt)
            if subarray and subarray != 'FULL':
                must_contain.append(subarray)
            else:
                for suba in MIRI_SUBARRAYS:
                    must_not_contain.append(suba)                      
            must_contain.append('DARK')
            if averaged:
                must_contain.append('averaged')
            else:
                must_not_contain.append('averaged')
            (matched_files, matched_folders) = \
                miricdp.match_cdp_substrings(mustcontain=must_contain,
                                             mustnotcontain=must_not_contain)
             
            if len(matched_files) > 1:
                # More than one match. Take the last file.
                filename = matched_files[-1]
                ftp_path = matched_folders[-1]
            elif len(matched_files) > 0:
                # Exactly one match. Take the first (and only) file.
                filename = matched_files[0]
                ftp_path = matched_folders[0]
            else:
                # No CDPs were matched
                filename = None
                 
            if filename:
                print("Matched %s at ftp_path=%s" % (filename, ftp_path))
                # Update the local CDP cache to make sure it contains the specified file,
                # and obtain the local file path and name.
                local_filename = miricdp.update_cache(filename, ftp_path)
                if averaged:
                    strg = "Reading averaged DARK model from \'%s\'" % local_filename
                else:
                    strg = "Reading DARK model from \'%s\'" % local_filename
                print(strg)
                dark_model = MiriDarkReferenceModel( init=local_filename )
            else:
                dark_model = None
            return dark_model

        print("Dark maps")
        for detector in ('MIRIMAGE', 'MIRIFUSHORT', 'MIRIFULONG'):
            for readpatt in ('FAST', 'SLOW'):
                for subarray in MIRI_SUBARRAYS:
                    strg = "DARK for detector=%s" % detector
                    strg += " readpatt=%s" % readpatt
                    strg += " subarray=%s" % subarray
                    print("\n" + strg)
                    darkmodel = find_dark( detector, readpatt, subarray,
                                           averaged=False)
                    if darkmodel is not None:
                        cdps_found.append(strg)
                        if VERBOSE_MODELS:
                            print( darkmodel )
                        if PLOTTING:
                            darkmodel.plot(strg)
                    else:
                        cdps_not_found.append(strg)
                        print("*** CDP NOT FOUND ***")
                    del darkmodel
                    time.sleep(LONG_DELAY)

        print("Pixel flat-fields")
        detector = 'MIRIMAGE'
        for detector in ('MIRIMAGE', 'MIRIFUSHORT', 'MIRIFULONG'):
            for mirifilter in ['ANY'] +  MIRI_FILTERS:
                for subarray in MIRI_SUBARRAYS:
                    strg = "PIXELFLAT for detector=%s" % detector
                    strg += " readpatt=%s" % readpatt
                    strg += " filter=%s" % mirifilter
                    strg += " subarray=%s" % subarray
                    print("\n" + strg)
                    flatmodel = get_cdp('PIXELFLAT', detector=detector,
                                        readpatt=readpatt, mirifilter=mirifilter,
                                        subarray=subarray,
                                        ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                        ftp_passwd=FTP_PASS)
                    if flatmodel is not None:
                        cdps_found.append(strg)
                        if VERBOSE_MODELS:
                            print( flatmodel )
                        if PLOTTING:
                            flatmodel.plot(strg)
                    else:
                        cdps_not_found.append(strg)
                        print("*** CDP NOT FOUND ***")
                    del flatmodel
                    time.sleep(LONG_DELAY)
 
        for detector in ('MIRIFUSHORT', 'MIRIFULONG'):
            for readpatt in ('ANY', 'FAST', 'SLOW'):
                strg = "PIXELFLAT for detector=%s" % detector
                strg += " readpatt=%s" % readpatt
                print("\n" + strg)
                flatmodel = get_cdp('PIXELFLAT', detector=detector,
                                    readpatt=readpatt,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if flatmodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( flatmodel )
                    if PLOTTING:
                        flatmodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del flatmodel
                time.sleep(LONG_DELAY)
 
        print("Sky flat-fields")
        detector = 'MIRIMAGE'
        for mirifilter in ['ANY'] +  MIRI_FILTERS:
            strg = "SKYFLAT for detector=%s" % detector
            strg += " filter=%s" % mirifilter
            print("\n" + strg)
            flatmodel = get_cdp('SKYFLAT', detector=detector,
                                mirifilter=mirifilter,
                                ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                ftp_passwd=FTP_PASS)
            if flatmodel is not None:
                cdps_found.append(strg)
                if VERBOSE_MODELS:
                    print( flatmodel )
                if PLOTTING:
                    flatmodel.plot(strg)
            else:
                cdps_not_found.append(strg)
                print("*** CDP NOT FOUND ***")
            del flatmodel
            time.sleep(LONG_DELAY)
 
        for detector in ('MIRIFUSHORT', 'MIRIFULONG'):
            for miriband in ['ANY'] + MIRI_BANDS:
                strg = "SKYFLAT for detector=%s" % detector
                strg += " band=%s" % miriband
                print("\n" + strg)
                flatmodel = get_cdp('SKYFLAT', detector=detector,
                                    band=miriband,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if flatmodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( flatmodel )
                    if PLOTTING:
                        flatmodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del flatmodel
                time.sleep(LONG_DELAY)
 
        print("Fringe flat-fields")
        detector = 'MIRIMAGE'
        for mirifilter in ['ANY'] +  MIRI_FILTERS:
            strg = "FRINGE for detector=%s" % detector
            strg += " filter=%s" % mirifilter
            print("\n" + strg)
            fringemodel = get_cdp('FRINGE', detector=detector,
                                mirifilter=mirifilter,
                                ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                ftp_passwd=FTP_PASS)
            if fringemodel is not None:
                cdps_found.append(strg)
                if VERBOSE_MODELS:
                    print( fringemodel )
                if PLOTTING:
                    fringemodel.plot(strg)
            else:
                cdps_not_found.append(strg)
                print("*** CDP NOT FOUND ***")
            del fringemodel
            time.sleep(LONG_DELAY)
 
        for detector in ('MIRIFUSHORT', 'MIRIFULONG'):
            for miriband in ['ANY'] + MIRI_BANDS:
                strg = "FRINGE for detector=%s" % detector
                strg += " band=%s" % miriband
                print("\n" + strg)
                fringemodel = get_cdp('FRINGE', detector=detector,
                                    band=miriband,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if fringemodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( fringemodel )
                    if PLOTTING:
                        fringemodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del fringemodel
                time.sleep(LONG_DELAY)
 
        print("Read noise models")
        detector = 'MIRIMAGE'
        for readpatt in ('FAST', 'SLOW'):
            for subarray in MIRI_SUBARRAYS:
                strg = "READNOISE for detector=%s" % detector
                strg += " readpatt=%s" % readpatt
                strg += " subarray=%s" % subarray
                print("\n" + strg)
                noisemodel = get_cdp('READNOISE', detector=detector,
                                    readpatt=readpatt, subarray=subarray,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if noisemodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( noisemodel )
                    if PLOTTING:
                        noisemodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del noisemodel
                time.sleep(LONG_DELAY)
 
        for detector in ('MIRIFUSHORT', 'MIRIFULONG'):
            for readpatt in ('FAST', 'SLOW'):
                strg = "READNOISE for detector=%s" % detector
                strg += " readpatt=%s" % readpatt
                print("\n" + strg)
                noisemodel = get_cdp('READNOISE', detector=detector,
                                    readpatt=readpatt,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if noisemodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( noisemodel )
                    if PLOTTING:
                        noisemodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del noisemodel
                time.sleep(LONG_DELAY)
 
        print("Distortion models")
        detector = 'MIRIMAGE'
        for mirifilter in ['ANY'] + MIRI_FILTERS:
            strg = "DISTORTION for detector=%s" % detector
            strg += " filter=%s" % mirifilter
            print("\n" + strg)
            distmodel = get_cdp('DISTORTION', detector=detector,
                                mirifilter=mirifilter,
                                ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                ftp_passwd=FTP_PASS)
            if distmodel is not None:
                cdps_found.append(strg)
                print( distmodel )
                if PLOTTING:
                    distmodel.plot(strg)
            else:
                cdps_not_found.append(strg)
                print("*** CDP NOT FOUND ***")
            del distmodel
            time.sleep(LONG_DELAY)
 
        for detector in ('MIRIFUSHORT', 'MIRIFULONG'):
            for miriband in ['ANY'] + MIRI_BANDS:
                strg = "DISTORTION for detector=%s" % detector
                strg += " band=%s" % miriband
                print("\n" + strg)
                distmodel = get_cdp('DISTORTION', detector=detector,
                                    band=miriband,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if distmodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( distmodel )
                    if PLOTTING:
                        distmodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del distmodel
                time.sleep(LONG_DELAY)
 
        print("Pixel area models")
        for detector in ('MIRIMAGE', 'MIRIFUSHORT', 'MIRIFULONG'):
            strg = "AREA for detector=%s" % detector
            print("\n" + strg)
            areamodel = get_cdp('AREA', detector=detector, ftp_path=FTP_PATH,
                                ftp_user=FTP_USER, ftp_passwd=FTP_PASS)
            if areamodel is not None:
                cdps_found.append(strg)
                if VERBOSE_MODELS:
                    print( areamodel )
                if PLOTTING:
                    areamodel.plot(strg)
            else:
                cdps_not_found.append(strg)
                print("*** CDP NOT FOUND ***")
            del areamodel
            time.sleep(LONG_DELAY)
 
        print("PSF models")
        detector = 'MIRIMAGE'
        for cdpmodel in ('PSF', 'PSF-OOF'):
            for mirifilter in ['ANY'] + MIRI_FILTERS:
                for subarray in MIRI_SUBARRAYS:
                    strg = "%s for detector=%s" % (cdpmodel, detector)
                    strg += " filter=%s" % mirifilter
                    strg += " subarray=%s" % subarray
                    print("\n" + strg)
                    psfmodel = get_cdp(cdpmodel, detector=detector,
                                        mirifilter=mirifilter, subarray=subarray,
                                        ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                        ftp_passwd=FTP_PASS)
                    if psfmodel is not None:
                        cdps_found.append(strg)
                        if VERBOSE_MODELS:
                            print( psfmodel )
                        if PLOTTING:
                            psfmodel.plot(strg)
                    else:
                        cdps_not_found.append(strg)
                        print("*** CDP NOT FOUND ***")
                    del psfmodel
                    time.sleep(LONG_DELAY)
         
        for detector in ('MIRIFUSHORT', 'MIRIFULONG'):
            for miriband in ('SHORT', 'MEDIUM', 'LONG'):
                strg = "PSF for detector=%s" % detector
                strg += " band=%s" % miriband
                print("\n" + strg)
                psfmodel = get_cdp('PSF', detector=detector,
                                    band=miriband,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if psfmodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( psfmodel )
                    if PLOTTING:
                        psfmodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del psfmodel
                time.sleep(LONG_DELAY)
 
        print("Photon Conversion Error (PCE) models")
        detector = 'MIRIMAGE'
        for mirifilter in ('ANY', 'P750L', 'F2100W'):
            for subarray in ('FULL', 'SLITLESSPRISM'):
                strg = "PCE for detector=%s" % detector
                strg += " filter=%s" % mirifilter
                strg += " subarray=%s" % subarray
                print("\n" + strg)
                pcemodel = get_cdp('PCE', detector=detector,
                                    mirifilter=mirifilter,
                                    subarray=subarray,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if pcemodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( pcemodel )
                    if PLOTTING:
                        pcemodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del pcemodel
                time.sleep(SHORT_DELAY)
 
        for detector in ('MIRIFUSHORT', 'MIRIFULONG'):
            for miriband in ('SHORT', 'MEDIUM', 'LONG'):
                strg = "PCE for detector=%s" % detector
                strg += " band=%s" % miriband
                print("\n" + strg)
                pcemodel = get_cdp('PCE', detector=detector,
                                    band=miriband,
                                    ftp_path=FTP_PATH, ftp_user=FTP_USER,
                                    ftp_passwd=FTP_PASS)
                if pcemodel is not None:
                    cdps_found.append(strg)
                    if VERBOSE_MODELS:
                        print( pcemodel )
                    if PLOTTING:
                        pcemodel.plot(strg)
                else:
                    cdps_not_found.append(strg)
                    print("*** CDP NOT FOUND ***")
                del pcemodel
                time.sleep(SHORT_DELAY)
                
        print("The following %d simulator CDPs were successfully found:" % \
              len(cdps_found))
        for found in cdps_found:
            print("\t" + found)
        time.sleep(SHORT_DELAY)
        print("The following %d simulator CDPs could not be found:" % \
              len(cdps_not_found))
        for notfound in cdps_not_found:
            print("\t" + notfound)
        time.sleep(SHORT_DELAY)
            
    if TEST_ERRORS:
        # The class must not crash if it can't connect to the CDP repository
        print("Test an authorisation error")
        miricdp = MiriCDPInterface(ftp_user='miri', ftp_passwd='xxxx')
        if VERBOSE_CDPIFACE:
            print( miricdp )
            # Wait for the print output catch up.
            time.sleep(MEDIUM_DELAY)
        del miricdp

        print("Test a file search error")
        miricdp = MiriCDPInterface(ftp_path='NOSUCHFOLDER', ftp_user=FTP_USER,
                                   ftp_passwd=FTP_PASS)
        if VERBOSE_CDPIFACE:
            print( miricdp )
            # Wait for the print output catch up.
            time.sleep(MEDIUM_DELAY)
        del miricdp
        
        print("Authorisation error while trying to get a CDP file")
        maskmodel = get_cdp('MASK', detector='MIRIMAGE', ftp_path='CDPSIM',
                            ftp_user=FTP_USER, ftp_passwd='xxxx')
        if VERBOSE_MODELS:
            print("Mask model returned=\n", maskmodel)
        time.sleep(MEDIUM_DELAY)
        del maskmodel

        print("File search error while trying to get a CDP file")
        maskmodel = get_cdp('MASK', detector='MIRIMAGE', ftp_path='NOSUCHFOLDER',
                            ftp_user=FTP_USER, ftp_passwd=FTP_PASS)
        if VERBOSE_MODELS:
            print("Mask model returned=\n", maskmodel)
        time.sleep(MEDIUM_DELAY)
        del maskmodel

        print("Attempt to get an unknown CDP file")
        puddingmodel = get_cdp('PUDDING', detector='MIRIMAGE', ftp_path='CDPSIM',
                            ftp_user=FTP_USER, ftp_passwd=FTP_PASS)
        if VERBOSE_MODELS:
            print("Pudding model returned=\n", puddingmodel)
        time.sleep(MEDIUM_DELAY)
        del puddingmodel

    print("Test finished.")
