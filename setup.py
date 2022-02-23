#!/usr/bin/env python

"""

Setup file for installing the MiriTE software

:History:
11 May 2017: Major new setup.py script based on ez_setup.py, using
             the mirisim setup.py as a template and copying code
             from the package defsetup.py scripts.
12 May 2017: MIRI software divided into 4 parts: MiriTools, MiriPipeline,
             MiriSimulators and MiriCalibration.
14 Jul 2017: Unpack the 470 micron version of the cosmic ray libraries.
27 Jul 2017: Require Python 2.7
13 Sep 2017: MiriCalibration and MiriPipeline packages separated.
15 Jan 2018: MiriTools and MiriSimulators levels removed from package.
22 Jan 2018: Removed empty mirimsim simulator package.
25 Jan 2018: Added missing url metadata.
27 Apr 2018: Require Python 3.5.
22 May 2018: Added README and LICENCE.
04 Jun 2018: Added warning if MIRICLE environment is not activated.
10 Aug 2018: Removed datamodel scripts which are no longer relevant.
08 Oct 2018: Added convert_mrs_resolution script.
11 Oct 2018: Added append_lrs_photom script.
12 Oct 2018: Added cdp_remove_junk script.
07 Oct 2019: Require Python 3.6. Corrected bug in the checking of
             conda_prefix.
23 Mar 2020: Require Python 3.7.
17 Apr 2020: MIRI-759: Restructured to make the developer install
             work the same as a regular install.
24 Apr 2020: Unzip the data files when creating a developer install.
             Stop unzipping obsolete detector files.
09 Jun 2020: Added make_sca_subarray.py
12 Jun 2020: Added 'install_requires' with required dependencies.
16 Jun 2020: Do not install obsolete amplifier calibration files.
             Do not install old and rarely used scrupts.
18 Jun 2020: Removed "numba" from dependencies due to issues with Numba,
             see MIRI-749.
01 Feb 2022: Restrict version of Cython dependency to avoid issue with conda.

@author: MIRI Software Team

"""

import io
import os
import re
import sys
import zipfile

try:
    from setuptools import Extension
#    from distutils.core import Extension
    from Cython.Distutils import build_ext
except ImportError:
    build_ext = None

try:
    from setuptools import setup
except ImportError:
    from .ez_setup import use_setuptools
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


def get_conda_prefix():
    import os
    if 'CONDA_PREFIX' in list(os.environ.keys()):
        conda_prefix = str(os.environ['CONDA_PREFIX'])
    else:
        conda_prefix = ''
    return conda_prefix


# Test the command arguments given with this script.
# Only unzip data files when building or installing, not when
# cleaning. If a "clean" is requested on its own, the previously
# unzipped files are deleted.
if len(sys.argv[0]) > 0:
    argv = sys.argv
if ("clean" in argv):
    zipflag = False
    cleanflag = not (("build" in argv) or ("install" in argv) or ("develop" in argv))
elif ("check" in argv):
    zipflag = False
    cleanflag = False
else:
    zipflag = True
    cleanflag = False
if "--quiet" in argv:
    verbose = False
else:
    verbose = True

conda_prefix = get_conda_prefix()
if verbose:
    print("CONDA_PREFIX is", conda_prefix)

# ------------------------------------------------------------------
# Unzip the data files contained in the simulators data directories.
#
# The simulator detector data files are found relative to the
# directory containing this Python script.
(this_dir, this_file) = os.path.split(__file__)
cr_data_path = os.path.join(this_dir, "miri/simulators/data/cosmic_rays")
if not os.path.isdir(cr_data_path):
    strg = "Cosmic ray data directory %s not found" % cr_data_path
    raise EnvironmentError(strg)

detector_data_path = os.path.join(this_dir, "miri/simulators/data/detector")
if not os.path.isdir(detector_data_path):
    strg = "Detector data directory %s not found" % detector_data_path
    raise EnvironmentError(strg)

scasim_data_path = os.path.join(this_dir, "miri/simulators/scasim/data")
if not os.path.isdir(scasim_data_path):
    strg = "SCASIM data directory %s not found" % scasim_data_path
    raise EnvironmentError(strg)

# Unpack the 470 micron version of the cosmic ray library (which has better track lengths)
fziplist0 = [zipfile.ZipFile(os.path.join(cr_data_path,'CRs_SiAs_470.zip'),'r')]
# fziplist0 = [zipfile.ZipFile(os.path.join(cr_data_path,'CRs_SiAs_35.zip'),'r')]
#fziplist1 = [zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsIM.zip'),'r'),
#             zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsLW.zip'),'r'),
#             zipfile.ZipFile(os.path.join(detector_data_path,'bad_pixelsSW.zip'),'r'),
#             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapIM.zip'),'r'),
#             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapLW.zip'),'r'),
#             zipfile.ZipFile(os.path.join(detector_data_path,'dark_mapSW.zip'),'r')]
fziplist1 = []
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


entry_points = dict(asdf_extensions=['miri_datamodel = miri.datamodels.miri_extension:MIRIExtension'])

# ------------------------------------------------------------------


setup(
    name="miri",
    version=find_version("miri/__init__.py"),
    description="MIRI tools, data models and simulator software",
    url="https://github.com/JWST-MIRI/MiriTE",
    author="MIRI European Consortium",
    author_email="mirisim@roe.ac.uk",
    license="See LICENCE file",
    platforms=["Linux", "Mac OS X"],
    python_requires='>=3.7',
    packages=['miri',
              'miri.tools', 'miri.tools.tests',
              'miri.datamodels', 'miri.datamodels.tests',
              'miri.simulators', 'miri.simulators.tests',
              'miri.simulators.scasim', 'miri.simulators.scasim.tests',
              'miri.apt_parser',
              ],
    package_dir={
                 'miri': 'miri',
                 'miri.tools': 'miri/tools/',
                 'miri.tools.tests': 'miri/tools/tests',
                 'miri.datamodels': 'miri/datamodels/',
                 'miri.datamodels.tests': 'miri/datamodels/tests',
                 'miri.simulators': 'miri/simulators/',
                 'miri.simulators.tests': 'miri/simulators/tests',
                 'miri.simulators.scasim': 'miri/simulators/scasim',
                 'miri.simulators.scasim.tests': 'miri/simulators/scasim/tests',
                 'miri.apt_parser':'miri/apt_parser',
                },
    package_data={'miri.tools': ['data/__init__.py'],
                  'miri.datamodels': ['schemas/*.yaml', 'data/*.fits',
                                      'data/*.txt', 'data/__init__.py'],
                  'miri.simulators': ['schemas/*.yaml', 'data/*.fits',
                                      'data/*.txt', 'data/__init__.py',
                                      'data/detector/*.fits',
                                      'data/detector/*.txt',
                                      'data/cosmic_rays/*.fits',
                                      'data/cosmic_rays/*.txt',
                                      'data/filters/*.fits',
                                      'data/filters/*.txt'],
                  'miri.simulators.scasim': ['data/SCATestInput80x64.fits',
                                             'data/SCATestHorseHead1024.fits',
                                             'data/__init__.py'],
                  },
    scripts=['miri_installation_check.py',
#             'miri/datamodels/scripts/append_lrs_photom.py',
             'miri/datamodels/scripts/cdp_add_filter_band.py',
             'miri/datamodels/scripts/cdp_add_history.py',
             'miri/datamodels/scripts/cdp_add_subarray.py',
#             'miri/datamodels/scripts/cdp_correct_band.py',
#             'miri/datamodels/scripts/cdp_correct_wildcard.py',
             'miri/datamodels/scripts/cdp_get_doc.py',
             'miri/datamodels/scripts/cdp_print.py',
             'miri/datamodels/scripts/cdp_reduce_dark.py',
             'miri/datamodels/scripts/cdp_remove_junk.py',
             'miri/datamodels/scripts/cdp_verify.py',
             'miri/datamodels/scripts/convert_fits_to_asdf.py',
#             'miri/datamodels/scripts/convert_mrs_resolution.py',
             'miri/datamodels/scripts/convert_slope_data.py',
             'miri/datamodels/scripts/find_me_another.py',
             'miri/datamodels/scripts/dqflags_examples.py',
             'miri/datamodels/scripts/multicdp_band.csh',
             'miri/datamodels/scripts/multicdp_filter_band.csh',
             'miri/datamodels/scripts/multicdp_remove_junk.csh',
             'miri/datamodels/scripts/multicdp_subarray.csh',
             'miri/datamodels/scripts/multicdp_verify.py',
#             'miri/datamodels/scripts/multicdp_wildcard.csh',
#             'miri/simulators/scasim/scripts/make_bad_pixel_mask.py',
#             'miri/simulators/scasim/scripts/make_fringe_map.py',
             'miri/simulators/scasim/scripts/make_sca_calibration.py',
             'miri/simulators/scasim/scripts/make_sca_file.py',
             'miri/simulators/scasim/scripts/make_sca_subarray.py',
             'miri/simulators/scasim/scripts/convert_exposure_data.py',
#             'miri/simulators/scasim/scripts/detector_latency_test.py',
             'miri/simulators/scasim/scripts/plot_exposure_data.py',
             'miri/simulators/scasim/scripts/scasim.py',
             ],
    data_files=[('', ['LICENCE', 'README'])],
    entry_points=entry_points,
    install_requires=[
        'Cython>=0.29.15,<=0.29.25',
        'jwst>=0.18.0',
        'matplotlib>=3.1.0',
        'numpy>=1.18.1',
        'parameterized>=0.7.0',
        'paramiko==2.6.0',
        'pysftp==0.2.9',
        'scipy>=1.4.1',
        'pytest>=6.2.0',
    ],
)

if not cleanflag:
    if not ('/envs' in conda_prefix):
        print("\n*** WARNING: MIRI software installed into the root environment! ***")
        print("If you didn't want to do this, remove the above package from site-packages, execute")
        print("\n\tsource activate <name-of-miricle-environment>")
        print("\nand try again.")
