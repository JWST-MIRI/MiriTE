#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines the MIRI
fringe frequencies model.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

20 Nov 2015: Created (using the transmission correction table as a template).
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
14 Nov 2018: Replaced 'ANY' with 'N/A'. Explicitly set table column units
             based on the tunit definitions in the schema.
15 Nov 2018: Removed redundant function.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""

import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
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
        self.meta.reftype = 'FRINGEFREQ'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type:
            self.meta.model_type = model_type        

        # This is a reference data model.
        self._reference_model()
        
        if fringefreq_table is not None:
            try:
                self.fringefreq_table = fringefreq_table
            except (ValueError, TypeError) as e:
                strg = "fringefreq_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
         
        # Copy the table column units from the schema, if defined.
        fringefreq_units = self.set_table_units('fringefreq_table')


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
        testfringefreq1.set_instrument_metadata(detector='MIRIFUSHORT',
                                         ccc_pos='OPEN', channel='N/A',
                                         band='N/A')
        testfringefreq1.set_subarray_metadata('FULL')
        testfringefreq1.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', useafter='2015-11-20',
                                           description='Test data')
        print(testfringefreq1)
        if PLOTTING:
            testfringefreq1.plot(description="testfringefreq1")
        if SAVE_FILES:
            testfringefreq1.save("test_fringefreq_model1.fits", overwrite=True)
        del testfringefreq1
        
    print("Test finished.")
