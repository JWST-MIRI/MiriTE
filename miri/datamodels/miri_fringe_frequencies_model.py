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
26 Mar 2020: Ensure the model_type remains as originally defined when saving
             to a file.
29 Apr 2021: Updated to work with the new CDP structure defined in the yaml file.
             Modified the form of the fringefreq_table parameter, now parsed and
             assigned to the relevant datamodel attributes, changed the inputs to
             the small test when run as main to work with new CDP class, updated
             docstring

@author: Steven Beard (UKATC), Vincent Geers (UKATC), Patrick Kavanagh (DIAS)

"""

import warnings
import numpy as np

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
        
    fringefreq_tables: list of list of tuples or numpy record array (optional)
        A list of 4 lists or recarrays, the first three containing the fringe
        frequencies of SHORT, MEDIUM, and LONG, and the fourth containing
        the maximum residual fringe amplitude array.
        For each, either: A list of tuples containing columns in the aperture
        fringe frequency or maximum amplitude tables;
        Or: A numpy record array containing the same information as above.
        The tables must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_fringe_frequencies_mrs.schema"
    fieldnames = ('Slice', 'ffreq', 'dffreq', 'min_nfringes', 'max_nfringes', 'min_snr', 'pgram_res')
    
    def __init__(self, init=None, fringefreq_table=None, **kwargs):
        """
        
        Initialises the MiriMrsFringeFrequenciesModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsFringeFrequenciesModel, self).__init__(init=init, **kwargs)

        # Data type is fringe frequencies.
        self.meta.reftype = 'FRINGEFREQ'
        # Initialise the model type
        self._init_data_type()       
        # This is a reference data model.
        self._reference_model()
        
        if fringefreq_table is not None:
            try:
                self.fringefreq_table_short = fringefreq_table[0]
                self.fringefreq_table_medium = fringefreq_table[1]
                self.fringefreq_table_long = fringefreq_table[2]
                self.max_amp = fringefreq_table[3]
            except (ValueError, TypeError) as e:
                strg = "fringefreq_table must be a list of 4 lists or recarrays (SHORT, MEDIUM, LONG, MAX_AMP)"
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
         
        # Copy the table column units from the schema, if defined.
        _ = self.set_table_units('fringefreq_table_short')

    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriMrsFringeFrequenciesModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriMrsFringeFrequenciesModel module.")
    
    PLOTTING = False
    SAVE_FILES = False

    fringefreqentry = [(101.0, np.array([2.9,0]), np.array([0.3,0]), np.array([5,0]),
                        np.array([30, 0]), np.array([10,0]), np.array([0.005,0])),
                       (102.0, np.array([2.9, 0]), np.array([0.3, 0]), np.array([5, 0]),
                        np.array([30, 0]), np.array([10, 0]), np.array([0.005, 0])),
                       (103.0, np.array([2.9, 0]), np.array([0.3, 0]), np.array([5, 0]),
                        np.array([30, 0]), np.array([10, 0]), np.array([0.005, 0])),
                       (104.0, np.array([2.9, 0]), np.array([0.3, 0]), np.array([5, 0]),
                        np.array([30, 0]), np.array([10, 0]), np.array([0.005, 0]))
                       ]

    # max_amp arrays
    wav = np.linspace(4.0, 10.0, 20)
    amp = np.ones(wav.shape) * 0.2
    maxampentry = list(zip(wav.tolist(), amp.tolist()))

    # test data
    fringefreqdata = [fringefreqentry, fringefreqentry, fringefreqentry, maxampentry]

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
