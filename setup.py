#!/usr/bin/env python

"""

Setup file for installing the /teams/miri software

:History:
11 May 2017: Major new setup.py script based on ez_setup.py, using
             the mirisim setup.py as a template and copying code
             from the package defsetup.py scripts.
12 May 2017: MIRI software divided into 4 parts: MiriTools, MiriPipeline,
             MiriSimulators and MiriCalibration.
14 Jul 2017: Unpack the 470 micron version of the cosmic ray libraries.
27 Jul 2017: Require Python 2.7
13 Sep 2017: MiriCalibration and MiriPipeline packages separated.

@author: Steven Beard (UKATC)

"""

import io
import os
import re
import sys
import zipfile
import numpy

try:
    from setuptools import Extension
#    from distutils.core import Extension
    from Cython.Distutils import build_ext
except ImportError:
    build_ext = None

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
cr_data_path = os.path.join(this_dir, "MiriSimulators/simulators/data/cosmic_rays")
if not os.path.isdir(cr_data_path):
    strg = "Cosmic ray data directory %s not found" % cr_data_path
    raise EnvironmentError(strg)

detector_data_path = os.path.join(this_dir, "MiriSimulators/simulators/data/detector")
if not os.path.isdir(detector_data_path):
    strg = "Detector data directory %s not found" % detector_data_path
    raise EnvironmentError(strg)

scasim_data_path = os.path.join(this_dir, "MiriSimulators/simulators/scasim/data")
if not os.path.isdir(scasim_data_path):
    strg = "SCASIM data directory %s not found" % scasim_data_path
    raise EnvironmentError(strg)

# Unpack the 470 micron version of the cosmic ray library (which has better track lengths)
fziplist0 = [zipfile.ZipFile(os.path.join(cr_data_path,'CRs_SiAs_470.zip'),'r')]
# fziplist0 = [zipfile.ZipFile(os.path.join(cr_data_path,'CRs_SiAs_35.zip'),'r')]
fziplist1 = [zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsIM.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsLW.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsSW.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapIM.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapLW.zip'),'r'),
             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapSW.zip'),'r')]
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


setup(
    name="miri",
    version=find_version("__init__.py"),
    description="MIRI Data Calibration and Analysis Software",
    author="MIRI European Consortium",
    author_email="F.Lahuis@sron.nl",
    license="See LICENCE file",
    platforms=["Linux", "Mac OS X", "Win"],
    python_requires='~=2.7',
    packages=['miri',
              'miri.tools', 'miri.tools.tests',
              'miri.datamodels', 'miri.datamodels.tests',
              'miri.simulators', 'miri.simulators.tests',
              'miri.simulators.mirimsim', 'miri.simulators.mirimsim.tests',
              'miri.simulators.scasim', 'miri.simulators.scasim.tests',
             ],
    package_dir={
                 'miri': '',
                 'miri.tools': 'MiriTools/tools/',
                 'miri.tools.tests': 'MiriTools/tools/tests',
                 'miri.datamodels': 'MiriTools/datamodels/',
                 'miri.datamodels.tests': 'MiriTools/datamodels/tests',
                 'miri.simulators': 'MiriSimulators/simulators/',
                 'miri.simulators.tests': 'MiriSimulators/simulators/tests',
                 'miri.simulators.mirimsim': 'MiriSimulators/simulators/mirimsim/',
                 'miri.simulators.mirimsim.tests': 'MiriSimulators/simulators/mirimsim/tests',
                 'miri.simulators.scasim': 'MiriSimulators/simulators/scasim',
                 'miri.simulators.scasim.tests': 'MiriSimulators/simulators/scasim/tests',
                },
    package_data={'miri.tools': ['data/__init__.py'],
                  'miri.datamodels': ['schemas/*.yaml', 'data/*.fits',
                                   'data/*.txt', 'data/__init__.py'],
                  'miri.simulators': ['schemas/*.yaml', 'data/*.fits',
                                      'data/*.txt', 'data/__init__.py',
                                      'data/amplifiers/*.fits',
                                      'data/amplifiers/*.txt',
                                      'data/detector/*.fits',
                                      'data/detector/*.txt',
                                      'data/cosmic_rays/*.fits',
                                      'data/cosmic_rays/*.txt',
                                      'data/filters/*.fits',
                                      'data/filters/*.txt'],
                  'miri.simulators.mirimsim': ['data/__init__.py'],
                  'miri.simulators.scasim': ['data/SCATestInput80x64.fits',
                                      'data/SCATestHorseHead1024.fits',
                                      'data/__init__.py'],
                 },
    scripts=['miri_installation_check.py',
             'MiriTools/datamodels/scripts/check_jwslib_datamodel.py',
             'MiriTools/datamodels/scripts/cdp_add_history.py',
             'MiriTools/datamodels/scripts/cdp_add_subarray.py',
             'MiriTools/datamodels/scripts/cdp_get_doc.py',
             'MiriTools/datamodels/scripts/cdp_print.py',
             'MiriTools/datamodels/scripts/cdp_reduce_dark.py',
             'MiriTools/datamodels/scripts/cdp_verify.py',
             'MiriTools/datamodels/scripts/convert_fits_to_asdf.py',
             'MiriTools/datamodels/scripts/find_me_another.py',
             'MiriTools/datamodels/scripts/dqflags_examples.py',
             'MiriTools/datamodels/scripts/multicdp_verify.py',
             'MiriTools/datamodels/scripts/multicdp_subarray.csh',
             'MiriTools/datamodels/scripts/make_filters_fits.py',
             'MiriTools/datamodels/scripts/make_measurements_fits.py',
             'MiriTools/datamodels/scripts/make_qe_fits.py',
             'MiriSimulators/simulators/mirimsim/scripts/mirimsim.py',
             'MiriSimulators/simulators/scasim/scripts/make_bad_pixel_mask.py',
             'MiriSimulators/simulators/scasim/scripts/make_fringe_map.py',
             'MiriSimulators/simulators/scasim/scripts/make_sca_calibration.py',
             'MiriSimulators/simulators/scasim/scripts/make_sca_file.py',
             'MiriSimulators/simulators/scasim/scripts/convert_exposure_data.py',
             'MiriSimulators/simulators/scasim/scripts/detector_latency_test.py',
             'MiriSimulators/simulators/scasim/scripts/plot_exposure_data.py',
             'MiriSimulators/simulators/scasim/scripts/scasim.py',
            ],
)
