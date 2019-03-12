#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module filesearching - Contains miscellaneous file searching functions.

Some of the file searching functions are based on ideas suggested by the
Python Cookbook, 2nd edition, by Alex Martelli, Anna Martelli Revenscroft
and David Ascher (O'Reilly Media 2005) 0_596-00797-3.

The class hierarchy is::

    ParameterFileManager

:History:

21 Oct 2011: Created
24 Oct 2011: Data directories included in default search path.
             add_to_searchpath and find_file_prefix functions added.
25 Oct 2011: ParameterFileManager class added.
31 Oct 2011: find_writable function added.
             Deprecated use of the .has_key() dictionary method removed.
11 Jan 2012: Renamed symbol to avoid potential name clash: file-->tryfile
16 Mar 2012: Corrected Sphinx documentation typos.
04 Apr 2012: Improvements suggested by pylint.
26 Apr 2012: Corrected the exception message raised by ParameterFileManager
             when the parameter file execution fails.
12 Nov 2012: _defaultpath modified so this module doesn't need a knowledge
             of the existence of other modules.
08 Sep 2015: Made compatible with Python 3. Changed default verbosity to
             reduce the amount of print output.
23 Mar 2016: Verbosity parameter replaced by Python logger.
30 Mar 2016: Added a search path parameter to the ParameterFileManager,
             so the search result can be made predictable.
26 Mar 2018: Changed 'nosuchfile' to a name even less likely to exist.
17 May 2018: Python 3: Converted dictionary keys return into a list.
:12 Mar 2019: Removed use of astropy.extern.six (since Python 2 no longer used).

@author: Steven Beard (UKATC)

"""

#from astropy.extern import six 
import sys, os, fnmatch

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.filesearch") # Get a default parent logger

# The default search path for MIRI utilities is the current directory
# followed by the data directories for all known MIRI modules.
import miri
_defaultpath = "." + os.pathsep + miri.__path__[0]


def add_to_searchpath( newpath, pathsep=os.pathsep ):
    """
    
    Add a specified directory to the global default file search path.
    
    NOTE: The global default file search path starts out as:
    
    * The current working directory.
    * The installed miri/data directory.
    * The installed miri.tools/data directory.
    * The installed miri/scasim/data directory.
    
    Further directories are added to the end of this path.
    Note also that this function affects the global default used when
    functions such as find_file_in_path or find_file are called without
    specifying a search path. It is possible to call these functions
    with any custom search path.
    
    :Parameters:
    
    newpath: str
        A string containing the name of the directory to be appended
        to the global default search path.
    pathsep: str, optional (default=os.pathsep)
        The string used to separate the directories in the search path.
        If not specified this defaults to os.pathsep, which is the path
        separator used by the current operating system (; for Windows
        and : for Linux).
    
    """
    global _defaultpath
    _defaultpath += pathsep + newpath

def make_searchpath( dir_list, pathsep=os.path.pathsep ):
    """
    
    Convert a list of directories into a search path string.
    
    
    :Parameters:
    
    dir_list: list of str
        A list of directories to convert into a search path.
    pathsep: str, optional (default=os.pathsep)
        The string used to separate the directories in the search path.
        If not specified this defaults to os.pathsep, which is the path
        separator used by the current operating system (; for Windows
        and : for Linux).
        
    :Returns:
    
    searchpath: str
        The search path. If the input list is empty an empty string is
        returned.
     
    
    """
    assert isinstance(dir_list, (tuple,list))
    if len(dir_list) < 1:
        return ''
    else:
        searchpath = pathsep.join( dir_list )
        return searchpath

def find_files_matching( root, patterns='*', single_level=False,
                         sortfiles=True, yield_folders=False,
                         yield_path_only=False):
    """
    
    Search all the subdirectories within a given root directory and
    generate a list of files whose names match the given pattern.

    Note: This is a generator rather than a function, so it may be
    used in Python loops like this::
    
        for file in find_files_matching('.', '*.fits'):
            print file
    
    A list may be obtained by calling this function within list().
    
    :Parameters:
    
    root: str
        The path of the root directory (e.g. '.' or '/tmp').
    patterns: str, optional
        A string containing a pattern to be matched. Multiple patterns
        may be given, separated by ';' (e.g. '*.fits', '*.py;*.htm;*.html').
        If no pattern is provided, all files will be matched.
    single_level: bool, optional (default=False)
        Set to True if only the top level directory is to be searched.
        If not specified, all nested subdirectories will be searched.
    sortfiles: bool, optional (default=True)
        Set to True if the list of files should be sorted before being
        yielded. This is the default. If set to False, the ordering
        returned by os.walk will be preserved.
    yield_folders: bool, optional (default=False)
        Set to True if folders matching the pattern are needed as well
        as simple files. If not specified, only files are yielded.
    yield_path_only: bool, optional (default=False)
        Set to True to yield the path to the folder containing the file
        rather than the file itself.
        
    :Yields:
    
    Next-filename: str
        The next matching file in the sequence
    
    """
    # Convert the semicolon separated patterns into a list.
    patterns = patterns.split(';')
    for path, subdirs, files in os.walk(root):
        # Add the folders to the list of files if needed.
        if yield_folders:
            files.extend(subdirs)
        # Optionally, sort the list of files before using it.
        if sortfiles:
            files.sort()
        # Scan through the list of files and list of patterns
        # and yield each file matching the pattern.
        for fname in files:
            for pattern in patterns:
                if fnmatch.fnmatch(fname, pattern):
                    if yield_path_only:
                        yield path
                    else:
                        yield os.path.join(path,fname)
                    break
        # If restricted to a single level, only the first path
        # is searched.
        if single_level:
            break

def find_file_in_tree( filename, root, path_only=False):
    """
    
    Search all the subdirectories within the given root directory
    and find the first file that matches the given file name.
    
    :Parameters:
    
    filename: str
        The name of the file to be located. The name can contain a
        pattern (e.g. 'miri*.txt').
    root: str
        The root directory from which to start searching for the file
        (e.g. '.' or '/tmp').
    path_only: bool, optional (default=False)
        Set to True to return the path to the folder containing the file
        rather than the path to the file itself.
        
    :Returns:
    
    filepath: str
        The full path (and name) of the matching file.
        If no file is found an empty string is returned.
    
    """
    for tryfile in find_files_matching(root, patterns=filename, sortfiles=False,
                                    yield_path_only=path_only):
        if tryfile:
            return os.path.abspath(tryfile)
    return ''

def _check_folder(folder, filename, walkdir=False, path_only=False):
    """
    
    Helper function to check for the presence of a file within
    a specified folder and optionally search for the file within
    that folder.
    
    """
    candidate = os.path.join(folder, filename)
    # The file will be recognised as a file if it exists.
    if os.path.isfile(candidate):
        if path_only:
            return os.path.abspath(folder)
        else:
            return os.path.abspath(candidate)
    # If the file is not found in the top level folder, the entire
    # directory tree underneath it can optionally be searched as well.
    if walkdir:
        matched = find_file_in_tree(filename, root=folder, path_only=path_only)
        if matched:
            return matched
    # No file found - return an empty string.
    return ''

def find_file_in_path( filename, search_path='.', pathsep=os.pathsep,
                       walkdir=False, path_only=False):
    """
    
    Find a file matching the given name in a given search path of
    directories.
    
    :Parameters:
    
    filename: str
        The name of the file to be located.
    search_path: str
        A search path string containing the directories to be searched.
        The string is assumed to contain a list of directories separated
        by a separator (e.g. '.:/bin:/local/bin')
    pathsep: str, optional (default=os.pathsep)
        The string used to separate the directories in the search path.
        If not specified this defaults to os.pathsep, which is the path
        separator used by the current operating system (; for Windows
        and : for Linux).
    walkdir: bool, optional (default=False)
        Set to True if the entire directory tree under each of the
        directories in the search path should be searched as well.
        By default, only the specific directories listed in the search
        path are searched.
        N.B. Specifying walkdir=True can slow down a file search
        significantly. It is best to specify only a short search path
        when walkdir=True.
    path_only: bool, optional (default=False)
        Set to True to return the path to the folder containing the file
        rather than the path to the file itself.
        
    :Returns:
    
    filepath: str
        The full path (and name) of the matching file.
        If no file is found an empty string is returned.
    
    """
    # Step through each path in the search path
    for path in search_path.split(pathsep):
        result = _check_folder(path, filename, walkdir=walkdir,
                               path_only=path_only)
        if result:
            return result
    # No file found - return an empty string.
    return ''
 
def find_file_in_pythonpath( filename, subfolder='miri',
                             walkdir=False, path_only=False):
    """
    
    Find a file matching the given name within the PYTHONPATH.
    
    :Parameters:
    
    filename: str
        The name of the file to be located.
    subfolder: str, optional (default='miri')
        If specified, a preferred subfolder within the PYTHONPATH which,
        if it exists, will be searched before the top level directory.
    walkdir: bool, optional (default=False)
        Set to True if the entire directory tree under each of the
        directories in the search path should be searched as well.
        By default, only the specific directories listed in the search
        path are searched.
        N.B. Specifying walkdir=True can slow down a file search
        significantly, especially if there are a lot of utilities
        installed in PYTHONPATH. The search can be speeded up by
        specifying a preferred subfolder.
    path_only: bool, optional (default=False)
        Set to True to return the path to the folder containing the file
        rather than the path to the file itself.
        
    :Returns:
    
    filepath: str
        The full path (and name) of the matching file.
        If no file is found an empty string is returned.
    
    """
    # Step through each path in the PYTHONPATH
    for path in sys.path:
        # If there is a preferred subfolder, and it exists
        # search this first.
        if subfolder:
            sfname = os.path.join(path, subfolder)
            if os.path.isdir(sfname):
                result = _check_folder(sfname, filename, walkdir=walkdir,
                                       path_only=path_only)
                if result:
                    return result
                
        result = _check_folder(path, filename, walkdir=walkdir,
                               path_only=path_only)
        if result:
            return result
    # No file found - return an empty string.
    return ''

def find_writable( search_path=None, pathsep=os.pathsep):
    """
    
    Find a writeable directory within the search path provided.
    By default, this function will search for a file within the current
    directory, the MIRI module data directories and the "/tmp" directory.
    
    :Parameters:
    
    search_path: str, optional (default=global default)
        A search path string containing the directories to be tested.
        The string is assumed to contain a list of directories separated
        by a separator (e.g. '.:/bin:/local/bin')
    pathsep: str, optional (default=os.pathsep)
        The string used to separate the directories in the search path.
        If not specified this defaults to os.pathsep, which is the path
        separator used by the current operating system (; for Windows
        and : for Linux).
        
    :Returns:
    
    path: str
        The first writable path found, or an empty string if none could
        be found.
    
    """
    global _defaultpath
    if search_path is None:
        search_path = _defaultpath + pathsep + "/tmp"
    for path in search_path.split(pathsep):
        if os.access(path, os.W_OK):
            return path
    # No writable directory found - return an empty string.
    return ''

def find_file( filename, search_path=None, pathsep=os.pathsep,
               walkdir=True, canraise=True):
    """
    
    Find a file matching the given name first in a given search path of
    directories and (failing that) within the PYTHONPATH.
    By default, this function will search for a file within the current
    directory the MIRI module data directories and then try PYTHONPATH.
    
    :Parameters:
    
    filename: str
        The name of the file to be located.
    search_path: str, optional (default=global default)
        A search path string containing the directories to be searched.
        The string is assumed to contain a list of directories separated
        by a separator (e.g. '.:/bin:/local/bin').
        By default, the current directory and data directories (plus
        any additional directories specified in add_to_search_path) are
        searched before moving on to PYTHONPATH.
    pathsep: str, optional (default=os.pathsep)
        The string used to separate the directories in the search path.
        If not specified this defaults to os.pathsep, which is the path
        separator used by the current operating system (; for Windows
        and : for Linux).
    walkdir: bool, optional (default=True)
        Set to True if the entire directory tree under each of the
        directories in the given search path should be searched as well.
        By default all directories are searched.
    canraise: bool, optional (default=True)
        Set to True if the function should fail with an exception if
        the given file is not found. Otherwise a failed search is
        indicated by the return of an empty string.
        
    :Returns:
    
    filepath: str
        The full path and name of the matching file.
        If no file is found an empty string is returned.    
    
    :Raises:
    
    NameError
        Raised if canraise=True and the file is not found.

    """
    global _defaultpath
    if search_path is None:
        search_path = _defaultpath
    
    firsttry = find_file_in_path( filename, search_path, pathsep=pathsep,
                                  walkdir=walkdir)
    if firsttry:
        return firsttry
    secondtry = find_file_in_pythonpath( filename, walkdir=walkdir )
    if secondtry:
        return secondtry
    else:
        if canraise:
            strg = "find_file: File \'%s\' not found in \n   \'%s\'" % \
                (filename, search_path)
            if walkdir:
                strg += " (in any subdirectory)"
            strg += " or in PYTHONPATH."
            raise NameError(strg)
        return ''

def find_file_prefix( fileprefix, search_path=None, pathsep=os.pathsep,
               walkdir=True, canraise=True):
    """
    
    Find the folder containing files with the given named prefix, first
    in a given search path of directories and (failing that) within the
    PYTHONPATH.
    By default, this function will search for a file within the current
    directory the MIRI module  data directories and then try PYTHONPATH.
    
    This function is used when, rather than there being a single parameter
    file, there are a whole set of parameter files with a common prefix.
    
    :Parameters:
    
    fileprefix: str
        A prefix for the files to be located (for example
        ''CR_SiAs_SUNMIN_'', which will cause all files matching
        ''CR_SiAs_SUNMIN_*.*'' to be located).
    search_path: str, optional (default=global default)
        A search path string containing the directories to be searched.
        The string is assumed to contain a list of directories separated
        by a separator (e.g. '.:/bin:/local/bin').
        By default, the current directory and data directories (plus
        any additional directories specified in add_to_search_path) are
        searched before moving on to PYTHONPATH.
    pathsep: str, optional (default=os.pathsep)
        The string used to separate the directories in the search path.
        If not specified this defaults to os.pathsep, which is the path
        separator used by the current operating system (; for Windows
        and : for Linux).
    walkdir: bool, optional (default=True)
        Set to True if the entire directory tree under each of the
        directories in the given search path should be searched as well.
        By default all directories are searched.
    canraise: bool, optional (default=True)
        Set to True if the function should fail with an exception if
        no matching files are found. Otherwise a failed search is
        indicated by the return of only the file prefix
        
    :Returns:
    
    fileprefixpath: str
        If matching files are found the full path of the folder containing
        the files plus the prefix is returned (e.g. ''/data/CR_SiAs_SUNMIN_'').
        If no files are found only the file prefix is returned (e.g.
        ''CR_SiAs_SUNMIN_'').
    
    :Raises:
    
    NameError
        Raised if canraise=True and no matching files are found.

    """
    global _defaultpath
    if search_path is None:
        search_path = _defaultpath

    filename = fileprefix + '*.*'
    firsttry = find_file_in_path( filename, search_path, pathsep=pathsep,
                                  walkdir=walkdir, path_only=True)
    if firsttry:
        return os.path.join(firsttry, fileprefix)
    secondtry = find_file_in_pythonpath( filename, walkdir=walkdir,
                                         path_only=True )
    if secondtry:
        return os.path.join(secondtry, fileprefix)
    else:
        if canraise:
            strg = "find_file: No files matching \'%s\'" % filename
            strg += " could be found in \n   \'%s\'" % search_path
            if walkdir:
                strg += " (in any subdirectory)"
            strg += " or in PYTHONPATH."
            raise NameError(strg)
        return fileprefix


class ParameterFileManager(object):
    """
    
    Class ParameterFileManager - searches for a named parameters file,
    processes it through the Python interpreter and stores all the
    variables defined within it. The class behaves rather like a Python
    dictionary, in that each Python variable defined in the file can be
    looked up by keyword. Multiple levels of keywords are supported, so
    the file may itself contain nested dictionaries and lists.
    
    :Parameters:
    
    filename: str
        The name of the parameters file to be located.     
    search_path: str
        A search path string containing the directories to be searched.
        The string is assumed to contain a list of directories separated
        by the OS path separator (e.g. '.:/bin:/local/bin').
        Defaults to the current directory.
    pathsep: str, optional (default=os.pathsep)
        The string used to separate the directories in the search path.
        If not specified this defaults to os.pathsep, which is the path
        separator used by the current operating system (; for Windows
        and : for Linux).
    description: str, optional, default="parameters"
        An optional description of the content of the parameters file.
    logger: Logger object (optional)
        A Python logger to handle the I/O. This parameter can be used
        by a caller to direct the output to a different logger, if
        the default defined by this module is not suitable.
    
    :Examples:
    
    Create a ParameterFileManager object from a particular file name::

        pfmobject = ParameterFileManager( 'parameter_file.py' )
    
    Get a list of the parameters defined within the parameter file.
    (Private and internal system variables are suppressed.) ::

        kwlist = pfmobject.keys()
        
    Check whether the ParameterFileManager object contains a keyword::

        result = keyword in pfmobject

    Get the value of a Python variable by specifying its name as a keyword::

        value = pfmobject[keyword]
        
    If the Python variable is itself a dictionary, look up a second
    keyword to get an entry from that dictionary. (This works with
    lists and tuples as well, except the second keyword should be an
    integer.) ::

        value = pfmobject.get(keyword1, keyword2)
    
    """
    # Store a list of the names of already opened files in a static cache.
    _filename_cache = []
        
    def __init__(self, filename, search_path='.', pathsep=os.pathsep,
                 description="parameters", logger=LOGGER):
        """
        
        Create a new ParameterFileManager object. 
                 
        Parameters: See class doc string.
       
        """
        self._logger = logger
        self._kwdict = { }
                  
        # First find the named parameter file somewhere within the
        # given search path. An empty string is returned if the file
        # is not found.
        fname = find_file_in_path(filename, search_path=search_path,
                                  pathsep=pathsep)
        if fname:
            # Report the full path and name of the file (important because
            # a user might have substituted their own file). However,
            # keep track of the file names in a cache so that each unique
            # file is reported only once.
            if not ParameterFileManager._filename_cache.__contains__(fname):
                ParameterFileManager._filename_cache.append(fname)
                self._logger.info( "Reading %s from parameter file %s" % \
                    (description,fname) )
            # Open the file, execute it with the Python interpreter and
            # store the result in self._kwdict.
            # NOTE: If there are compilation errors (as might happen if a
            # user introduces a typo into a parameter file) an exception
            # will be raised here. The file must always be closed.
            fp = open(fname,"r")
            try:
                code = fp.read()
                #six.exec_(code, self._kwdict)
                exec(code, self._kwdict)
            except Exception as e:
                strg = "Error while creating ParameterFileManager object.\n"
                strg += "  Failed to interpret %s as a Python file.\n" % fname
                strg += "  %s: %s" % (e.__class__.__name__, e)
                raise Exception(strg)
            finally:
                try:
                    fp.close()
                except Exception:
                    pass
            self._fname = fname
        else:
            strg = "Error while creating ParameterFileManager object.\n"
            strg += "  Parameter file %s not found." % filename
            raise NameError(strg)
        
    def has_key(self, keyword):
        """
        
        Returns True if the given keyword is contained in the
        ParameterFileManager object.
                 
        :Parameters:
    
        keyword: str
            The keyword to be tested.
        
        :Returns:
        
        has_key: bool
            True if the keyword exists, False if it does not.
         
        """
        #NOTE dictionary.has_key() is deprecated.
        return (keyword in self._kwdict)

    def keys(self):
        """
        
        Returns a list of keywords known to the ParameterFileManager
        object (ignoring private keywords beginning with '__' or system
        keywords associated with modules, functions and classes).
                 
        :Parameters:
    
        None.
        
        :Returns:
        
        keylist: list of str
            List of known keywords
        
        """
        keylist = []
        for key in list(self._kwdict.keys()):
            # Ignore keywords beginning '__' and module and function
            # definitions containing '<module', '<function' or '<class.
            if key.find('__') != 0:
                value = str(self._kwdict[key])
                if value.find('<module') == -1 and \
                   value.find('<function') == -1 and \
                   value.find('<class') == -1:
                    keylist.append(key)
        return keylist
    
    def __contains__(self, keyword):
        """
        
        Returns True if the given keyword is known to the
        ParameterFileManager object.
        
        """
        return (keyword in self._kwdict)
        
    def __setitem__(self, keyword, value):
        """
        
        Set an entry with the operation::

            pfmobject[keyword] = newvalue
            
        THIS IS NOT ALLOWED.
        
        """
        strg = "Modifying keyword %s in ParameterFileManager" % keyword
        strg += " is not allowed."
        raise TypeError(strg)
    
    def __delitem__(self, keyword):
        """
        
        Remove the given keyword entry from the ParameterFileManager object
        with the operation::

            del pfmobject[keyword]
            
        THIS IS NOT ALLOWED.
        
        """
        strg = "Removal of keyword %s from ParameterFileManager" % keyword
        strg += " is not allowed."
        raise TypeError(strg)
         
    def get(self, keyword1, keyword2=None, keyword3=None):
        """
        
        Return the value associated with the given set of keywords.
        Up to 3 levels of keywords can be specified, which allows
        values contained in nested dictionaries, tuples and lists to
        be looked up in one function call.
                 
        :Parameters:
    
        keyword1: str
            The keyword with which to look up a top level variable
            stored in the ParameterFileManager object.
        keyword2: str or int, optional (default=None)
            If the item looked up by keyword1 is a nested dictionary,
            tuple or list, subitems can be lookup up directly by
            specifying a second keyword.
        keyword3: str or int, optional (default=None)
            If the item looked up by keyword2 is a nested dictionary,
            tuple or list, subitems can be lookup up directly by
            specifying a third keyword.
 
        :Returns:

        value: python object
            The value of the variable defined in the parameter file.
    
        """
        #TODO: Nested keywords could be handled recursively.
        level2 = self._kwdict[keyword1]
        if keyword2 is None or not keyword2:
            return level2
        else:
            if isinstance(level2, dict):
                level3 = level2[str(keyword2)]
            elif isinstance(level2, (tuple,list)):
                level3 = level2[int(keyword2)]
            else:
                strg = "Item %s is not a dictionary, tuple or list." % keyword1
                strg += " Cannot look up subitem with %s." % keyword2
                raise TypeError(strg)
            
            if keyword3 is None or not keyword3:
                return level3
            else:
                if isinstance(level3, dict):
                    level4 = level3[str(keyword3)]
                elif isinstance(level3, (tuple,list)):
                    level4 = level3[int(keyword3)]
                else:
                    strg = "Item %s.%s is not a dictionary, tuple or list." % \
                        (keyword1, keyword2)
                    strg += " Cannot look up subitem with %s." % keyword3
                    raise TypeError(strg)
                return level4
   
    def __getitem__(self, keyword):
        """
        
        Return the value associated with the given keyword with
        the operation::

            value = pfmobject[keyword]
                    
        """
        return self._kwdict[keyword]
    
    def __str__(self):
        """
        
        Return a string description of the parameters contained
        within a ParameterFileManager object.
        
        """
        strg = ''
        strg += "ParameterFileManager: \n"
        strg += "---------------------\n"
        strg += "Source file: %s\n" % self._fname
        # List the values assiciated with each of the known
        # non-private keywords. If a keyword is another dictionary,
        # the values contained inside that are listed as well.
        for key in list(self.keys()):
            value = self._kwdict[key]
            if isinstance(value, dict):
                strg += "%8s :\n" % key
                for key2 in list(value.keys()):
                    value2 = value[key2]
                    strg += "\t%8s (" % key2
                    strg += "\'%-8s\')\n" % str(value2)
            else:
                strg += "%8s (" % key
                strg += "\'%-8s\')\n" % str(value)
        return strg


if __name__ == '__main__':
    print( "Testing the filesearching module." )

    test_files = ('filesearching.py', \
                  'miriplot.py', \
                  '__init__.py', \
                  'example_measurement.txt', \
                  'bad_pixelsLW.fits', \
                  'nequetamlimanoexisteelarchivo42.txt')

    for name in test_files:
        fmatched = find_file(name, canraise=False)
        if fmatched:
            print( "%s found at %s" % (name, fmatched) )
        else:
            print( "%s not found" % name )

    name = 'CRs_SiAs_SUNMIN_'
    fmatched = find_file_prefix(name, canraise=False)
    if fmatched != name:
        print( "%s found at %s" % (name, fmatched) )
    else:
        print( "%s not found" % name )

    print( "" )
    # BEWARE: example_properties.py imports this module. The only reason
    # there isn't an infinite recursion loop is because the import skips
    # this __main__ part.
    import miri.datamodels
    dir_list = ['.',
                os.path.dirname(__file__),
                miri.datamodels.__path__[0]]
    search_path = make_searchpath(dir_list)     
    example_properties = ParameterFileManager("example_properties.py",
                                        search_path=search_path,
                                        description="some example properties")
    print( example_properties )

    print( "Test finished." )
