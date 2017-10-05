#!/usr/bin/env python

"""

Setup file for installing the MIRI simulator software

:History:
:History:
02 Jun 2017: New stand-alone setup script for the MIRI simulator
             software only.

@author: Steven Beard (UKATC)

"""

import io
import os
import re
import sys
import zipfile
import numpy

try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

    
def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


# Test the command arguments given with this script.
# Only unzip data files when building or installing, not when 
# cleaning. If a "clean" is requested on its own, the previously 
# unzipped files are deleted.
if len(sys.argv[0]) > 0:
    argv = sys.argv
if ("build" in argv) or ("install" in argv):
    zipflag = True
    cleanflag = False
else:
    zipflag = False
    cleanflag = ("clean" in argv)
if "--quiet" in argv:
    verbose = False
else:
    verbose = True


# ------------------------------------------------------------------
# Unzip the data files contained in the simulators data directories.
#
# The simulator detector data files are found relative to the 
# directory containing this Python script.
(this_dir, this_file) = os.path.split(__file__)

cr_data_path = os.path.join(this_dir, "simulators/data/cosmic_rays")
if not os.path.isdir(cr_data_path):
    strg = "Cosmic ray data directory %s not found" % cr_data_path
    raise EnvironmentError(strg)

detector_data_path = os.path.join(this_dir, "simulators/data/detector")
if not os.path.isdir(detector_data_path):
    strg = "Detector data directory %s not found" % detector_data_path
    raise EnvironmentError(strg)

scasim_data_path = os.path.join(this_dir, "simulators/scasim/data")
if not os.path.isdir(scasim_data_path):
    strg = "SCASIM data directory %s not found" % scasim_data_path
    raise EnvironmentError(strg)

fziplist0 = [zipfile.ZipFile(os.path.join(cr_data_path,'CRs_SiAs_35.zip'),'r')]
fziplist1 = [zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsIM.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsLW.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsSW.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapIM.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapLW.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapSW.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'MIRI_FM_MIRIMAGE_FAST_DARK.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'MIRI_FM_MIRIFULONG_FAST_DARK.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'MIRI_FM_MIRIFUSHORT_FAST_DARK.zip'),'r')]
fziplist2 = [zipfile.ZipFile(os.path.join(scasim_data_path,'SCATestHorseHead1024.zip'),'r')]

if zipflag:
    # Unpack the cosmic ray files
    for fzip in fziplist0:
        for name in fzip.namelist():
            fullname = os.path.join(cr_data_path, name)
            if not os.path.isfile(fullname):
                if verbose:
                    print( "Unzipping \'%s\'" % fullname )
                data = fzip.read(name)
                temp = open(fullname, "wb")
                temp.write(data)
                temp.close()
            else:
                if verbose:
                    print( "\'%s\' already exists" % fullname )
        fzip.close()
    # Unpack the detector files
    for fzip in fziplist1:
        for name in fzip.namelist():
            fullname = os.path.join(detector_data_path, name)
            if not os.path.isfile(fullname):
                if verbose:
                    print( "Unzipping \'%s\'" % fullname )
                data = fzip.read(name)
                temp = open(fullname, "wb")
                temp.write(data)
                temp.close()
            else:
                if verbose:
                    print( "\'%s\' already exists" % fullname )
        fzip.close()
    # Unpack the SCASIM files
    for fzip in fziplist2:
        for name in fzip.namelist():
            fullname = os.path.join(scasim_data_path, name)
            if not os.path.isfile(fullname):
                if verbose:
                    print( "Unzipping \'%s\'" % fullname )
                data = fzip.read(name)
                temp = open(fullname, "wb")
                temp.write(data)
                temp.close()
            else:
                if verbose:
                    print( "\'%s\' already exists" % fullname )
        fzip.close()

elif cleanflag:
    # Clean up the cosmic ray files
    for fzip in fziplist0:
        for name in fzip.namelist():
            fullname = os.path.join(cr_data_path, name)
            if os.path.isfile(fullname):
                if verbose:
                    print( "Deleting \'%s\'" % fullname )
                try:
                    os.remove(fullname)
                except Exception:
                    pass
    # Clean up the detector files
    for fzip in fziplist1:
        for name in fzip.namelist():
            fullname = os.path.join(detector_data_path, name)
            if os.path.isfile(fullname):
                if verbose:
                    print( "Deleting \'%s\'" % fullname )
                try:
                    os.remove(fullname)
                except Exception:
                    pass
    # Clean up the SCASIM files
    for fzip in fziplist2:
        for name in fzip.namelist():
            fullname = os.path.join(scasim_data_path, name)
            if os.path.isfile(fullname):
                if verbose:
                    print( "Deleting \'%s\'" % fullname )
                try:
                    os.remove(fullname)
                except Exception:
                    pass

# ------------------------------------------------------------------

# TODO: Change the scasim.simulators.scasim package structure to just scasim
setup(
    name="scasim",
    version=find_version("__init__.py"),
    description="MIRI Simulator Software",
    author="MIRI European Consortium",
    author_email="steven.beard@stfc.ac.uk",
    license="See LICENCE file",
    platforms=["Linux", "Mac OS X", "Win"],
    python_requires='~=2.7',
    packages=['scasim',
              'scasim.simulators', 'scasim.simulators.tests',
              'scasim.simulators.mirimsim', 'scasim.simulators.mirimsim.tests',
              'scasim.simulators.scasim', 'scasim.simulators.scasim.tests',
             ],
    package_dir={
                 'scasim': '',
                 'scasim.simulators': 'simulators/',
                 'scasim.simulators.tests': 'simulators/tests',
                 'scasim.simulators.mirimsim': 'simulators/mirimsim/',
                 'scasim.simulators.mirimsim.tests': 'simulators/mirimsim/tests',
                 'scasim.simulators.scasim': 'simulators/scasim',
                 'scasim.simulators.scasim.tests': 'simulators/scasim/tests',
                },
    package_data={'scasim.simulators': ['schemas/*.yaml', 'data/*.fits',
                                      'data/*.txt', 'data/__init__.py',
                                      'data/amplifiers/*.fits',
                                      'data/amplifiers/*.txt',
                                      'data/detector/*.fits',
                                      'data/detector/*.txt',
                                      'data/cosmic_rays/*.fits',
                                      'data/cosmic_rays/*.txt',
                                      'data/filters/*.fits',
                                      'data/filters/*.txt'],
                  'scasim.simulators.mirimsim': ['data/__init__.py'],
                  'scasim.simulators.scasim': ['data/SCATestInput80x64.fits',
                                      'data/SCATestHorseHead1024.fits',
                                      'data/__init__.py'],
                 },
    scripts=['simulators/mirimsim/scripts/mirimsim.py',
             'simulators/scasim/scripts/make_bad_pixel_mask.py',
             'simulators/scasim/scripts/make_fringe_map.py',
             'simulators/scasim/scripts/make_sca_calibration.py',
             'simulators/scasim/scripts/make_sca_file.py',
             'simulators/scasim/scripts/convert_exposure_data.py',
             'simulators/scasim/scripts/detector_latency_test.py',
             'simulators/scasim/scripts/plot_exposure_data.py',
             'simulators/scasim/scripts/scasim.py',
            ],
)
