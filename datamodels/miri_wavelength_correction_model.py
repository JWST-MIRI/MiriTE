#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines the MIRI
wavelength correction model.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

30 Jun 2015: Created
08 Jul 2015: Changed column headings in optical table.
09 Jul 2015: Removed duplication of table units between schema and metadata.
             Units are now only defined in the metadata.
             Use of the fieldnames class variable removed from the code and
             deprecated. It is now used only by a few conversion scripts.
11 Sep 2015: Removed duplicated plot method.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""
# This module is now converted to Python 3.


#import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriMrsWavelengthCorrectionModel']


class MiriMrsWavelengthCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI wavelength correction table,
    based on the STScI base model, DataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape.
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    wavcorr_optical: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the wavelength
        correction optical table;
        Or: A numpy record array containing the same information as above.
        A wavelength correction table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    wavcorr_xslice: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the wavelength
        correction cross slice table;
        Or: A numpy record array containing the same information as above.
        A wavelength correction table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    wavcorr_shift: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the wavelength
        correction spectral shift table;
        Or: A numpy record array containing the same information as above.
        A wavelength correction table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_wavelength_correction_mrs.schema.yaml"
    fieldnames_optical = ('SUB_BAND', 'BETA_SLICE', 'WAVE_MIN', 'WAVE_MAX',
                          'SRP_MIN', 'SRP_MAX')
    fieldnames_xslice =  ('XSLICE_MIN', 'XSLICE_MAX')
    fieldnames_shift =   ('BETA_OFF', 'DS_MIN', 'DS_MAX')
    
    def __init__(self, init=None, wavcorr_optical=None, wavcorr_xslice=None,
                 wavcorr_shift=None, **kwargs):
        """
        
        Initialises the MiriMrsWavelengthCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsWavelengthCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is wavelength correction.
        self.meta.model_type = 'WAVCORR'
        self.meta.reftype = 'WAVCORR'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'

        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'
        
        if wavcorr_optical is not None:
            try:
                self.wavcorr_optical = wavcorr_optical
            except (ValueError, TypeError) as e:
                strg = "wavcorr_optical must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if wavcorr_xslice is not None:
            try:
                self.wavcorr_xslice = wavcorr_xslice
            except (ValueError, TypeError) as e:
                strg = "wavcorr_xslice must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if wavcorr_shift is not None:
            try:
                self.wavcorr_shift = wavcorr_shift
            except (ValueError, TypeError) as e:
                strg = "wavcorr_shift must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        
#         # Copy the table column units, if defined.
#         wavcorr_optical_units = self.set_table_units('wavcorr_optical')
#         wavcorr_xslice_units = self.set_table_units('wavcorr_xslice')
#         wavcorr_shift_units = self.set_table_units('wavcorr_shift')
        
    def __str__(self):
        """
        
        Return the contents of the wavelength correction object
        as a readable string.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')

        # Describe the wavelength correction tables
        if self.wavcorr_optical is not None:
            strg += self.get_data_str('wavcorr_optical', underline=True, underchar="-")
        if self.wavcorr_xslice is not None:
            strg += self.get_data_str('wavcorr_xslice', underline=True, underchar="-")
        if self.wavcorr_shift is not None:
            strg += self.get_data_str('wavcorr_shift', underline=True, underchar="-")
        return strg


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriMrsWavelengthCorrectionModel module.")
    
    PLOTTING = False
    SAVE_FILES = False

    optical_data = [('1A', 0.176, 4.87, 5.82, 3320.0, 3710.0),
                    ('1B', 0.176, 5.62, 6.73, 3190.0, 3750.0)]
    xslice_data =  [(0.2122, 0.3718)]
    shift_data =   [(0.000, 0.0, 0.0),
                    (0.005, -0.0460, -0.0687),
                    (0.010, -0.0924, -0.0687)]
    
    print("\nWavelength correction with factors derived from list of tuples:")
    with MiriMrsWavelengthCorrectionModel( wavcorr_optical=optical_data,
                                           wavcorr_xslice=xslice_data,
                                           wavcorr_shift=shift_data ) \
                                           as testwavcorr1:
        print(testwavcorr1)
        if PLOTTING:
            testwavcorr1.plot(description="testwavcorr1")
        if SAVE_FILES:
            testwavcorr1.save("test_wavcorr_model1.fits", overwrite=True)
        del testwavcorr1
        
    print("Test finished.")
