#!/usr/bin/env python

"""

Setup file for installing the MIRI data models and software tools.

:History:
15 May 2017: New script which builds only the MiriTools part of the
             MIRI software.
02 Jun 2017: Corrected the details.

@author: Steven Beard (UKATC)

"""

import io
import os
import re
import sys

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


setup(
    name="miri",
    version=find_version("__init__.py"),
    description="MIRI Data Models and Software Tools",
    author="MIRI European Consortium",
    author_email="steven.beard@stfc.ac.uk",
    license="See LICENCE file",
    platforms=["Linux", "Mac OS X", "Win"],
    python_requires='~=2.7',
    packages=['miri',
              'miri.tools', 'miri.tools.tests',
              'miri.datamodels', 'miri.datamodels.tests',
             ],
    package_dir={
                 'miri': '',
                 'miri.tools': 'tools/',
                 'miri.tools.tests': 'tools/tests',
                 'miri.datamodels': 'datamodels/',
                 'miri.datamodels.tests': 'datamodels/tests',
                },
    package_data={'miri.tools': ['data/__init__.py'],
                  'miri.datamodels': ['schemas/*.yaml', 'data/*.fits',
                                   'data/*.txt', 'data/__init__.py'],
                 },
    scripts=['datamodels/scripts/check_jwslib_datamodel.py',
             'datamodels/scripts/cdp_add_history.py',
             'datamodels/scripts/cdp_add_subarray.py',
             'datamodels/scripts/cdp_get_doc.py',
             'datamodels/scripts/cdp_print.py',
             'datamodels/scripts/cdp_reduce_dark.py',
             'datamodels/scripts/cdp_verify.py',
             'datamodels/scripts/convert_fits_to_asdf.py',
             'datamodels/scripts/find_me_another.py',
             'datamodels/scripts/dqflags_examples.py',
             'datamodels/scripts/multicdp_verify.py',
             'datamodels/scripts/multicdp_subarray.csh',
             'datamodels/scripts/make_filters_fits.py',
             'datamodels/scripts/make_measurements_fits.py',
             'datamodels/scripts/make_qe_fits.py',
            ],
)
