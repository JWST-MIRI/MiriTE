#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines the MIRI
fringe frequencies model.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

20 Nov 2015: Created (using the transmission correction table as a template).
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""
# This module is now converted to Python 3.


import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriMrsFringeFrequenciesModel']


class MiriMrsFringeFrequenciesModel(MiriDataModel):
    """
    
    A generic data model for a MIRI fringe frequencies table,
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
        
    fringefreq_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the aperture
        correction table;
        Or: A numpy record array containing the same information as above.
        A fringe frequencies table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_fringe_frequencies_mrs.schema.yaml"
    fieldnames = ('sub_band', 'wavenumber', 'deltawavenumber', 'maxamplitude')
    
    def __init__(self, init=None, fringefreq_table=None, **kwargs):
        """
        
        Initialises the MiriMrsFringeFrequenciesModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsFringeFrequenciesModel, self).__init__(init=init, **kwargs)

        # Data type is fringe frequencies.
        self.meta.model_type = 'FRINGEFREQ'
        self.meta.reftype = 'FRINGEFREQ'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'
        
        if fringefreq_table is not None:
            try:
                self.fringefreq_table = fringefreq_table
            except (ValueError, TypeError) as e:
                strg = "fringefreq_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
#         
#         # Copy the table column units, if defined.
#         fringefreq_units = self.set_table_units('fringefreq_table')
        
    # TODO: Is this function needed?
    def __str__(self):
        """
        
        Return the contents of the fringe frequencies object
        as a readable string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()

        # Describe the fringe frequencies table
        if self.fringefreq_table is not None:
            strg += self.get_data_str('fringefreq_table', underline=True, underchar="-")
        return strg


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriMrsFringeFrequenciesModel module.")
    
    PLOTTING = False
    SAVE_FILES = False

    fringefreqdata = [('ANY', 11.0,  0.1, 42.0),
                      ('ANY', 12.0,  0.2, 22.0),
                      ('ANY', 13.0,  0.3, 32.0),
                      ('ANY', 14.0,  0.4, 12.0)]

    print("\nFringe frequencies with factors derived from list of tuples:")
    with MiriMrsFringeFrequenciesModel( fringefreq_table=fringefreqdata ) as testfringefreq1:
        print(testfringefreq1)
        testfringefreq1.set_instrument_metadata(detector='MIRIFUSHORT',
                                         ccc_pos='OPEN', channel='ANY',
                                         band='ANY')
        testfringefreq1.set_subarray_metadata('FULL')
        testfringefreq1.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', useafter='2015-11-20',
                                           description='Test data')
        if PLOTTING:
            testfringefreq1.plot(description="testfringefreq1")
        if SAVE_FILES:
            testfringefreq1.save("test_fringefreq_model1.fits", overwrite=True)
        del testfringefreq1
        
    print("Test finished.")
