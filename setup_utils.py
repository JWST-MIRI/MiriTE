#-*- coding: utf-8 -*-
#!/usr/bin/env python

"""

Package setup_utils provides utilities for setup.py

== OBSOLETE MODULE (replaced by ez_setup.py) ==

:History:

13 Dec 2010: Created with functions checkDeps, cleanInstall, rmtree.
31 Oct 2011: Modified rmtree so it can cope with symbolic links.
12 Mar 2012: Corrected typos.
04 Apr 2012: Some tidying up.
10 Apr 2012: Added the checkPythonpath utility.
12 Apr 2012: Use the clean and install options together to remove a
             previous installation instead of asking the question.
06 Nov 2012: Moved some inline code from setup.py into here.
             Protect file reads against exceptions with try/finally
             blocks.
07 Feb 2013: Updated list of dependencies. Warnings issued using the
             Python warnings module (so the messages are logged rather
             than lost when when the screen scrolls up).
28 Nov 2013: exec_package_defsetups and write_defsetup added here to
             simplify the top level setup.py.
02 Dec 2013: Avoid using the reserved word "license" and the Python
             deprecated term "licence". Added some "assert" statements
             to catch faulty user input.
06 Dec 2013: Catch other exceptions besides ImportError when checking
             dependent modules. Simplified the PYTHONPATH check.
08 Sep 2015: Made compatible with Python 3 (except for unicode literals,
             which cause a problem in the interpretation of package names).
07 Oct 2015: This module is not compatible with Python 3! See FIXME.
20 Sep 2016: Removed dependency on astropy.
12 May 2017: all_packages_input now gives the name and location of each
             package.

@author: Julien Morin (DIAS), Steven Beard (UKATC)

"""
from __future__ import absolute_import, division, print_function

#from io import open

import os
import sys
import warnings
try:
    from astropy.extern import six
    string_types = six.string_types
except ImportError:
    six = None
    string_types = (str, unicode)


def check_pythonpath( givenpath=sys.path ):
    """
    
    Checks the integrity of the given search path (default sys.path).
    It checks that all folders mentioned in the path actually exist
    and looks for duplicated entries (which can sometimes lead to import
    problems).
    
    """
    strg = ''
    # First check all the references in the PYTHONPATH actually exist.
    missing = []
    for pfolder in givenpath:
        if not os.path.exists(pfolder):
            missing.append(pfolder)
    if missing:
        strg += "\nPYTHONPATH contains the following non-existent entries:\n"
        strg += "  " + "  ".join(missing)
    
    # Then check for duplicated entries.
    duplicated = []
    for pfolder in givenpath:
        if givenpath.count(pfolder) > 1:
            duplicated.append(pfolder)      
    if duplicated:
        strg += "\nPYTHONPATH contains the following duplicate entries:\n"
        strg += "  " + "  ".join(duplicated)
    # Combine both messages into one warning.
    if strg:
        strg += "\nCorrect these PYTHONPATH issues if you see unexpected " \
            "ImportError or AccessInit errors."
        warnings.warn(strg)

def check_deps(deps, optdeps):
    """

    Check if the python packages listed in deps (compulsory dependencies)
    and optdeps (optional dependencies) are present in the PYTHONPATH.
    If compulsory dependencies are missing an exception is raised.
    If optional dependencies are missing a warning is issued.

    """
    # Both deps and optdeps must be a list of strings (or empty lists)
    assert isinstance(deps, (tuple,list))
    if deps:
        assert isinstance(deps[0], string_types)
    assert isinstance(optdeps, (tuple,list))
    if deps:
        assert isinstance(optdeps[0], string_types)

    deps_notfound = []
    optdeps_notfound = []
    for lib in deps:
        try:
            exec("import %s" % str(lib))
        except ImportError:
            deps_notfound.append(str(lib))
        except Exception as e:
            strg = "Unexpected exception while importing module %s:" % str(lib)
            strg += "\n  %s. " % str(e)
            strg += "\n  Installation aborted."
            raise e.__class__(strg)
    for lib in optdeps:
        try:
            exec("import %s" % str(lib))
        except ImportError:
            optdeps_notfound.append(str(lib))
        except Exception as e:
            strg = "Unexpected exception while importing module %s:" % str(lib)
            strg += "\n  %s: %s. " % (e.__class__.__name__, str(e))
            strg += "\n  Some functionality may not be available."
            warnings.warn(strg)

    if len(deps_notfound): 
        strg = "The following compulsory dependent packages are missing:\n"
        strg += "  " + ", ".join(deps_notfound)
        strg += "\n  Installation aborted. "
        strg += "Please install these packages and try again."
        raise ImportError(strg)
    elif len(optdeps_notfound): 
        strg = "The following optional dependent packages are missing:\n"
        strg += "  " + ", ".join(optdeps_notfound)
        strg += "\n  Some functionality may not be available."
        warnings.warn(strg)

def clean_install( ):
    """

    Check for previous miri install present in the PYTHONPATH. If a
    miri package is found it is removed. Although not strictly required
    this procedure ensures that unnecessary files are removed and *.pyc
    files updated.
    
    """
    # NOTE: Using force=1 in setup.cfg might make this function unnecessary.
    old_miri_path = None
    try :
        import miri as old_miri
        old_miri_path = old_miri.__path__[0]
        print( "  Previous miri install found at %s" % old_miri_path )
    except (ImportError, NameError, IndexError) :
        print( "  No previous miri install was found." )

    if old_miri_path:
        rmtree(old_miri_path)
        print( "  %s cleaned." % old_miri_path )

def rmtree(path):
    """

    Recursively remove the content of a given path.

    NOTE: shutil.rmtree provides a similar functionality but does not work
    correctly on all platforms and Python releases.

    """
    if os.path.isdir(path):
        files_list = os.listdir(path)
        for targetfile in files_list:
            file_path = os.path.join(path, targetfile)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                rmtree(file_path)
                os.rmdir(file_path)
    else:
        return

def exec_defsetup( filename, argv=None ):
    """

    Parse and execute the contents of the given defsetup.py file and
    return the setupargs and pkg variables declared within it.
    
    The optional argv variable can be used to pass command line
    arguments to the defsetup.py script (which allows it to distinguish
    between a "clean" and an "install").

    """
    # Initialise the return variables.
    pkg = None
    setupargs = None

    # "import" the defsetup file.  We can't use import directly
    # because none of the defsetup.py files are in modules, but we
    # can read the file and exec its contents.
    # This assumes we are running the build from the top level miri
    # directory.
    #
    # We want to protect our current environment from contamination
    # by the exec'ed code, so we give it a private symbol table

    if os.path.isfile(filename) :
        e_dict = { }
        e_dict['argv'] = argv # Pass over the system arguments
        e_dict['__file__'] =  filename # Pass over the file name.

        try:
            fp = open(filename, "r")
            code = fp.read()
            if six:
                six.exec_(code, e_dict)
            else:
                exec(code, e_dict)
        finally:
            fp.close()
    else :
        strg = "defsetup.py file, %s, cannot be found." % filename
        strg += "\n    Installation aborted."
        raise Exception(strg)

    # Pick out the two interesting variables from the exec'ed code
    # and return them.
    if "pkg" in e_dict :
        pkg = e_dict["pkg"]
    if "setupargs" in e_dict :
        setupargs = e_dict["setupargs"]

    if isinstance(pkg, str) :
        pkg = [ pkg ]

    return (setupargs, pkg)

def exec_package_defsetups( namespace, all_packages_input, argv=None):
    """
    
    Parse all the defsetup.py files belonging to subpackages
    and return a combined list of package directories, scripts
    and data files, etc...

    The optional argv variable can be used to pass command line
    arguments to each defsetup.py script.
    
    Returns (all_packages, all_package_dir, all_ext_modules,
             all_scripts, all_data_files)
    
    """
    assert isinstance(namespace, string_types)
    assert isinstance(all_packages_input, dict)
    # Start the list of all packages found.  all_packages_input only lists 
    # things at the top level, but if one of them has nested packages, 
    # all_packages will be the complete list. Begin with the top level 
    # package.
    all_packages = [namespace]

    # Start the dictionary of package directories with the top level directory.
    all_package_dir = {namespace:''}

    # Start the with empty lists for native C code modules, scripts and
    # data files.
    all_ext_modules = [ ]
    all_scripts = [ ]
    all_data_files = [ ]

    # For each package we want to use, fetch the config data from defsetup.py.
    # We have to fix various things because each package description thinks 
    # it is only describing itself, but now we are in the parent directory 
    # and need all the references prefixed with the parent package name.
    #
    for lpkg in all_packages_input.keys() :
        #print("Package", lpkg, "location", all_packages_input[lpkg])

        # Begin by appending the namespace reference for this package to
        # the all_packages list.
        nlpkg = namespace + '.' + lpkg
        all_packages.append(nlpkg)

        # Execute the defsetup.py file for this package and extract the
        # setupargs and pkg variables defined within it.
        fname = os.path.join(all_packages_input[lpkg], "defsetup.py")
        print( "Executing", all_packages_input[lpkg], "package defsetup.py file:", fname )
        (setupargs, pkg) = exec_defsetup(fname, argv=argv)

        # If the package doesn't report the same name that we asked
        # for, there is a major problem.
        if nlpkg != pkg[0] :
            strg = "Package name expected (%s) and " % nlpkg
            strg += "in defsetup.py (%s) do not match." % pkg[0]
            strg += "\n    Installation aborted."
            raise Exception(strg)    

        # Append the package directories defined in defsetup.py to the 
        # all_package_dir list. Each directory needs to be prefixed with 
        # lpkg to make it referable from the top level namespace.
        if not 'package_dir' in setupargs :
            # If package_dir is not specified, use the top-level as default
            all_package_dir[lpkg] = os.path.join( all_packages_input[lpkg], '' )

        else :
            # There is a list of package dirs to handle
            package_dir = setupargs['package_dir']
            for x in package_dir :
                if not x in all_packages :
                    all_packages.append(x)
                all_package_dir[x] = os.path.join( all_packages_input[lpkg], package_dir[x] )

        # If there are scripts, we have to correct the file names where
        # the installer can find them.  Each is under the package directory.
        if 'scripts' in setupargs :
            for x in setupargs['scripts'] :
                all_scripts.append(os.path.join( all_packages_input[lpkg], x ))

        # If there are external modules, we need to correct the file names
        # of source files.
        if 'ext_modules' in setupargs :
            for x in setupargs['ext_modules'] :
                x.sources = prefix_dir_list( x.sources, all_packages_input[lpkg] )
                all_ext_modules.append(x)

        # If there are data files, we need to correct the file names.
        if 'data_files' in setupargs :
            for x in setupargs['data_files'] :
                ( instdir, files ) = x
                t = prefix_dir_list( files, all_packages_input[lpkg] )
                all_data_files.append( ( instdir, t ) )
                
    return (all_packages, all_package_dir, all_ext_modules, \
            all_scripts, all_data_files)

def write_defsetup( namespace, version, description, author, email, licence,
                    all_packages, all_package_dir,
                    all_ext_modules, all_scripts, all_data_files ):
    """
    
    Write a new defsetup.py file containing the information contained in
    the function arguments.
    
    NOTE: Any existing defsetup.py file will be overwritten.
    
    """
    df = open("defsetup.py", "w")
    try:
        df.write(str('pkgname = \"%s\"\n\n' % namespace))
        df.write(str('pkg = %s\n' % str(all_packages)))

        setupargs = '{\n'
        setupargs += '\t\"version\" :          \"%s",\n' % version
        setupargs += '\t\"description\" :      \"%s\",\n' % description
        setupargs += '\t\"author\" :           \"%s\",\n' % author
        setupargs += '\t\"maintainer_email\" : \"%s\",\n' % email
        setupargs += '\t\"license\" :          \"%s\",\n' % licence
        setupargs += '\t\"platforms\" :        [\"Linux\", \"Solaris\", \"Mac OS X\", \"Win\"],\n'
        setupargs += '\t\"package_dir\" :      %s,\n' % str(all_package_dir)
        setupargs += '\t\"ext_modules\" :      %s,\n' % str(all_ext_modules)
        setupargs += '\t\"scripts\" :          %s,\n' % str(all_scripts)
        setupargs += '\t\"data_files\" :       %s,\n' % str(all_data_files)
        setupargs += '}'
        
        df.write('setupargs = %s\n' % str(setupargs))

    finally:
        df.close()

def prefix_dir_list( inlist, prefix, separator=os.path.sep ):
    """

    Add the given prefix onto each entry in the given directory
    list and return the new list.

    """
    outlist = []
    for item in inlist :
        if not prefix or item.startswith( separator ) :
            # No prefix given, or item begins with separator.
            outlist.append( item )
        else :
            # Add the prefix.
            outlist.append( os.path.join( prefix, item ))
    return outlist

def file_copy( infile, outfile ):
    """

    Make a copy of the given file.

    """
    try:
        fp = open( infile, 'r' )
        content = fp.readlines()
    finally:
        fp.close()
    try:
        fp = open( outfile, 'w' )
        fp.writelines( content )
    finally:
        fp.close()
    del content

def file_append(filename, addedlines):
    """

    Append extra lines to the end of the given file

    """
    try:
        fp = open( filename, 'a' )
        fp.write( addedlines )
    finally:
        fp.close()

def file_add_package_imports(filename, all_packages_input):
    """

    Append package imports to the given Python file.

    """
    add_lines = ""
    for pkg in all_packages_input :
        add_lines += "try:\n"
        add_lines += "    import %s\n" % pkg
        add_lines += "    PACKAGES.append('%s')\n" % pkg
        add_lines += "except ImportError:\n"
        add_lines += "    print('Package %s could not be loaded')\n\n" % pkg
    file_append( filename, add_lines )

def file_script_list(filename, scriptlist, separator=os.path.sep):
    """

    Write a list of scripts to the given file.

    """
    try:
        fp = open(filename,"w")
        for x in scriptlist :
            x = x.split( separator )
            fp.write(x[0] + " " + x[-1] + "\n")
    finally:
        fp.close()

def file_collection(filename, collection, separator=os.path.sep):
    """

    Write a list or dictionary to the given file.

    """
    try:
        fp = open(filename,"w")
        if isinstance(collection, dict ) :
            keys = collection.keys()
            for key in keys :
                fp.write(key + " = " + str(collection[key]) + "\n")
        elif isinstance(collection, (tuple,list) ) :
            for x in collection :
                fp.write(str(x) + "\n")
    finally:
        fp.close()
