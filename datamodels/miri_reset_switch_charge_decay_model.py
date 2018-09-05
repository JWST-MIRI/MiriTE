#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An experimental extension to the standard STScI data model, which 
defines MIRI reset switch charge decay model.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

19 Nov 2015: Created as a skeleton. Three pieces of information needed:
             (1) Names of decay constants; (2) units of decay constants;
             and (3) any additional metadata.
27 Nov 2015: Updated to match actual data file delivered for CDP-5.
01 Jun 2016: "ROWS" column added to RSCD table.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


#import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriResetSwitchChargeDecayModel']


class MiriResetSwitchChargeDecayModel(MiriDataModel):
    """
    
    A generic data model for a MIRI RSCD table, based on the STScI
    base model, DataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the flux parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    rscd_table. list of tuples or numpy record array (optional)
        A list of tuples, or a numpy record array, where each record
        contains the following information:
        
        * SUBARRAY: A string containing a subarray name. The string can
          be 'GENERIC' to indicate that the decay factors do not
          depend on the subarray.
        * READPATT: A string containing a readout pattern. The string can
          be 'ANY' to indicate that the decay factors do not
          depend on the readout pattern.
        * ONE: First decay constant (TBD).
        * TWO: First decay constant (TBD).
        * THREE: First decay constant (TBD).
        * FOUR: First decay constant (TBD).

        A rscd_table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
     detector: str (optional)
        The name of the detector associated with this RSCD data.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_reset_switch_charge_decay.schema.yaml"
    fieldnames = ('subarray', 'readpatt', 'rowtype', 'tau', 'ascale', 'pow',
                  'illum_zp', 'illum_slope', 'illum2', 'param3',
                  'crossopt', 'sat_zp', 'sat_slope', 'sat_2', 'sat_mzp',
                  'sat_rowterm', 'sat_scale')
   
    def __init__(self, init=None, rscd_table=None, detector=None, **kwargs):
        """
        
        Initialises the MiriResetSwitchChargeDecayModel class.
        
        Parameters: See class doc string.

        """
        super(MiriResetSwitchChargeDecayModel, self).__init__(init=init, **kwargs)

        # Data type is RSCD.
        self.meta.model_type = 'RSCD'
        self.meta.reftype = 'RSCD'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'
        
        # Define the detector identifier, if specified. (N.B. The detector
        # ID is compulsory, so it must be specified either here, in the
        # source file or later after creation of this data object.)
        # The warning is commented out because it causes unnecessary
        # output during tests or when creating an empty data product and
        # filling it in later.
        if detector is not None:
            self.meta.instrument.detector = detector

        if rscd_table is not None:
            try:
                self.rscd_table = rscd_table
            except (ValueError, TypeError) as e:
                strg = "rscd_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
#         
#         # Copy the table column units, if defined.
#         rscd_units = self.set_table_units('rscd_table')
        
    # TODO: Is this function needed?
    def __str__(self):
        """
        
        Return the contents of the RSCD object as a readable
        string.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title_and_metadata()

        # Describe the RSCD table
        if self.rscd_table is not None:
            strg += self.get_data_str('rscd_table', underline=True, underchar="-")
        else:
            strg += "No RSCD table!"
        return strg


class MiriResetSwitchChargeDecayModel_CDP6(MiriResetSwitchChargeDecayModel):
    """
    
    This class can be used to access the old CDP-6 version of the
    MiriMrsDistortionModel12 data model.
    
    See the MiriMrsDistortionModel12 class for full documentation.
    
    """
    schema_url = "miri_reset_switch_charge_decay_CDP6.schema.yaml"
    fieldnames = ('SUBARRAY', 'READPATT', 'ROWS', 'TAU1', 'SCALE1', 'TAU2', 'SCALE2')
    
    def __init__(self, init=None, rscd_table=None, detector=None, **kwargs):
        """
        
        Initialises the MiriResetSwitchChargeDecayModel class.
        
        Parameters: See class doc string.

        """
        super(MiriResetSwitchChargeDecayModel_CDP6, self).__init__(init=init,
                                rscd_table=rscd_table, detector=detector,
                                **kwargs)
 

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriResetSwitchChargeDecayModel module.")
    
    PLOTTING = False
    SAVE_FILES = False
     
    rscddata = [('FULL',          'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('FULL',          'FAST', 'EVEN', 1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('FULL',          'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('FULL',          'SLOW', 'EVEN', 1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('MASK1065',      'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('MASK1065',      'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('MASK1140',      'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('MASK1140',      'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('MASK1550',      'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('MASK1550',      'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('MASKLYOT',      'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('MASKLYOT',      'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('BRIGHTSKY',     'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('BRIGHTSKY',     'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('SLITLESSPRISM', 'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ('SLITLESSPRISM', 'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4, 5.0e-5, 6.0e-6, 7.0e-7, 8.0e-8, 9.0e-9, 10.0e-2, 11.0e-4, 12.0e-6, 13.0e-8, 14.0),
                ]

    print("\nRSCD calibration with factors derived from list of tuples:")
    with MiriResetSwitchChargeDecayModel( rscd_table=rscddata,
                                          detector='MIRIFUSHORT' ) as testrscd1:
        testrscd1.set_instrument_metadata(detector='MIRIFUSHORT',
                                          ccc_pos='OPEN',
                                          deck_temperature=11.0,
                                          detector_temperature=6.0)
        testrscd1.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', useafter='2015-11-20',
                                           description='Test data')
        print(testrscd1)
        if PLOTTING:
            testrscd1.plot(description="testrscd1")
        if SAVE_FILES:
            testrscd1.save("test_rscd_model1.fits", overwrite=True)
        del testrscd1
        
    print("Test finished.")
