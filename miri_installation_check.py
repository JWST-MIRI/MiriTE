#!/usr/bin/env python

"""

This module checks a MIRI Python software installation, makes sure all
important modules can be imported and reports the versions and locations
of those modules (when such information is available).

02 Apr 2012: Created.
11 Apr 2012: check_module function implemented.
25 Jun 2013: Added astronomical Python modules.
04 Dec 2013: Removed dependence on PACKAGES variable.
11 Feb 2014: Check for jwst_lib and print warnings using warnings module.
24 Feb 2014: Renamed to miri_installation_check.py so it may be safely
             installed as a script.
27 May 2015: Replaced pyfits with astropy.io.fits
08 Sep 2015: Made compatible with Python 3.
16 Sep 2016: Removed unused sandbox and stpipeline packages.
27 Sep 2016: Added json and yaml (and asdf commented out).
28 Sep 2016: Renamed miri.miritools to miri.tools and miri.miritools.datamodels
             to miri.datamodels. Uncommented asdf check.
02 Jun 2017: Modified to match latest dependencies.
04 Sep 2017: Added memory test
13 Sep 2017: Removed MiriCalibration and MiriPipeline.
05 Jan 2018: SVN info dropped.

@author: Steven Beard (UKATC)

"""
from __future__ import absolute_import, unicode_literals, division, print_function

import warnings

# Minimum recommended memory for using MIRI data models. Most data models
# will work when less memory is available, but a workstation with less than
# this limit will struggle to process 4-D DARK and ramp data.
MBYTES_MIN = 8000.0

def check_module( mymodule, indent='' ):
    """
    
    Returns a string describing the version number, import path and
    other administrative information contained in a given module.
    
    :Parameters:
    
    mymodule: Python module
        The module to be checked. It must already have been imported.
        
    indent: str, optional, default=''
        A string used to pad the left side of the report for indentation
        
    :Returns:
    
    moduleinfo: str
        A string describing the module.
    
    """
    _strg = "%s%s" % (indent, repr(mymodule))
    if hasattr(mymodule, '__name__'):
        _strg += "\n%s   - name:  %s" % (indent, str(mymodule.__name__))
    if hasattr(mymodule, '__project__'):
        _strg += "\n%s   - project:  %s" % (indent, str(mymodule.__project__))
    if hasattr(mymodule, '__package__'):
        _strg += "\n%s   - package:  %s" % (indent, str(mymodule.__package__))
    if hasattr(mymodule, '__author__'):
        _strg += "\n%s   - author:  %s" % (indent, str(mymodule.__author__))
    if hasattr(mymodule, '__version__'):
        _strg += "\n%s   - version:  %s" % (indent, str(mymodule.__version__))
    if hasattr(mymodule, '__path__'):
        _strg += "\n%s   - imported from:  %s" % \
            (indent, str(mymodule.__path__))
    if hasattr(mymodule, '__file__'):
        _strg += "\n%s   - initialised from:  %s" % \
            (indent, str(mymodule.__file__))
    if hasattr(mymodule, '__all__'):
        nall = len(mymodule.__all__)
        _strg += "\n%s   - length of _all_ list: %d" % (indent, nall)
    if hasattr(mymodule, '__doc__'):
        if mymodule.__doc__ is not None:
            _strg += "\n%s   - (%d character doc string)" % \
                (indent, len(mymodule.__doc__))
        else:
            _strg += "\n%s   - (null doc string)" % indent
    else:
        _strg += "\n%s   - (no doc string)" % indent

    return _strg


if __name__ == '__main__':
    print( "Installation check starting" )
    print( "===========================" )
    installation_ok = True

    # Display information about the Python environment
    print( "\nPlatform and environment" )
    print( "------------------------" )
    import sys, os
    print( "sys.platform =", sys.platform )
    print( "sys.version =", sys.version )
    print( "sys.api_version =", sys.api_version )
    print( "sys.hexversion = %x" % sys.hexversion )

    # Display the total memory available and warn if smaller than recommended.
    memok = True
    try:
        import psutil
        result = psutil.virtual_memory()
        mbytes = result[0] / (1024.0 * 1024.0)
        strg = "Total memory available = %.1f MB" % mbytes
        if mbytes < MBYTES_MIN:
           memok = False
           strg += " (*Less than recommended minimum of %.0f MB!*)" % MBYTES_MIN
        strg += "\n"
        print( strg )
    except ImportError:
        pass

    # Display the PYTHONPATH
    strg = "PYTHONPATH ="
    for pfolder in sys.path:
        strg += "\n  %s%s" % (pfolder, os.pathsep)
    print( strg )

    # Attempt to import various important modules and display version numbers
    # and other administrative information
    # TODO: The same exception checks are repeated over and over again. Reuse?
    # It isn't possible to put these in a function because the module has to
    # be explicitly imported here.
    print( "\nPython modules" )
    print( "----------------" )
    try:
        import math, random, re, glob, datetime
    except ImportError as e:
        installation_ok = False
        print( "*** Basic modules math/random/re/glob/datetime could not be imported: %s" % e )
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing basic modules math/random/re/glob/datetime:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)
    try:
        import yaml
        print( check_module( yaml ) )
    except ImportError as e:
        installation_ok = False
        print( "*** Module yaml could not be imported: %s" % e )
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing module yaml:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)
    try:
        import paramiko
        print( check_module( paramiko ) )
        import pysftp
        print( check_module( pysftp ) )
    except ImportError as e:
        installation_ok = False
        print( "*** Modules paramiko/pysftp could not be imported: %s" % e )
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing modules paramiko/pysftp:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)

    print( "\nScientific Python modules" )
    print( "-------------------------" )
    try:
        import numpy
        print( check_module( numpy ) )
    except ImportError as e:
        installation_ok = False
        print( "*** Module numpy could not be imported: %s" % e )
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing module numpy:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)

    print( " " )
    try:
        import scipy
        print( check_module( scipy ) )
    except ImportError as e:
        installation_ok = False
        warnings.warn("*** Module scipy could not be imported: %s" %e)
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing module scipy:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)

    print( " " )
    try:
        import matplotlib
        strg = check_module( matplotlib )
        import matplotlib.pyplot as plt
        strg += "\n" + check_module( plt, indent='   ' )
        import matplotlib.image as mpimg
        strg += "\n" + check_module( mpimg, indent='   ' )
        print( strg )
    except ImportError as e:
        installation_ok = False
        warnings.warn("*** module matplotlib could not be imported: %s" % e)
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing module matplotlib:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)

    print( "\nAstronomical Python modules" )
    print( "-----------------------------" )
    try:
        import astropy
        strg = check_module( astropy )
        # The following checks depend on having numpy V1.7 or later
        import astropy.extern as astroextern
        strg += "\n" + check_module( astroextern, indent='   ' )
        import astropy.io.ascii as asciitable
        strg += "\n" + check_module( asciitable, indent='   ' )
        import astropy.io.fits as astrofits
        strg += "\n" + check_module( astrofits, indent='   ' )
        import astropy.wcs as astrowcs
        strg += "\n" + check_module( astrowcs, indent='   ' )
        print( strg )
    except ImportError as e:
        installation_ok = False
        warnings.warn("*** Module astropy could not be imported: %s" % e)
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing module astropy:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)

    print( "\nSTScI Python modules" )
    print( "--------------------" )
    try:
        import asdf
        print( check_module( asdf ) )
    except ImportError as e:
        installation_ok = False
        warnings.warn("*** Module asdf could not be imported: %s" % e)
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing module asdf:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)

    print( " " )
    try:
        import jwst.datamodels
        print( check_module( jwst.datamodels ) )
    except ImportError as e:
        installation_ok = False
        print( "*** Module jwst.datamodels could not be imported: %s" % e )
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing module jwst.datamodels:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)

    # Finally, check that the MIRI software can be imported.
    print( "\nMIRI software modules" )
    print( "---------------------" )
    try:
        import miri
        print( check_module( miri ) )
        print( " " )
        import miri.tools
        print( check_module( miri.tools, indent='   ') )
        print( " " )
        import miri.datamodels
        print( check_module( miri.datamodels, indent='   ') )
        print( " " )
        import miri.simulators
        print( check_module( miri.simulators, indent='   ') )
    except ImportError as e:
        installation_ok = False
        warnings.warn("*** Module miri could not be imported: %s" % e)
    except Exception as e:
        installation_ok = False
        strg = "*** Miscellaneous error while importing module miri:"
        strg += "\n  %s: %s" % (e.__class__.__name__, e)
        warnings.warn(strg)

    strg = "\nInstallation check finished. "
    if installation_ok and memok:
        strg += "Installation OK."
        print( strg )
    elif not installation_ok:
        strg += "Installation NOT OK - see problems reported above."
        warnings.warn(strg)
    elif not memok:
        strg += "Installation OK but available memory less than %.0f MB." % MBYTES_MIN
        warnings.warn(strg)
