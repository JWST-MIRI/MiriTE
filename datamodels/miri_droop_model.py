#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An experimental extension to the standard STScI data model, which 
defines MIRI droop model.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

22 Feb 2013: Created
17 May 2013: Do not allow a blank table to be created.
22 Aug 2013: columnnames renamed to fieldnames.
02 Sep 2013: Pass the responsibility for creating record arrays to jwst.datamodels
             - a solution to the "Types in column 0 do not match" problem
             suggested by Michael Droettboom at STScI.
03 Sep 2013: DETECTOR field removed from table in anticipation of CDP-2
             delivery. Check there is a DETECTOR identification in the
             header.
12 Sep 2013: Make the creation of an empty table a warning rather than an
             error.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: TYPE and REFTYPE are no longer identical.
09 Dec 2014: Removed old commented out warnings.
09 Jul 2015: Removed duplication of table units between schema and metadata.
             Units are now only defined in the metadata.
             Use of the fieldnames class variable removed from the code and
             deprecated. It is now used only by a few conversion scripts.
11 Sep 2015: Removed duplicated plot method.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

#import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriDroopModel']


class MiriDroopModel(MiriDataModel):
    """
    
    A generic data model for a MIRI droop table, based on the STScI
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
        
    droop_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (coupling_constant:number,
        uncertainty:number), giving the droop couping constant.
        Or: A numpy record array containing the same information as above.
        A droop table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
    detector: str (optional)
        The name of the detector associated with this droop data.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_droop.schema.yaml"
    fieldnames = ('COUPLING_CONSTANT', 'UNCERTAINTY')
    
    def __init__(self, init=None, droop_table=None, detector=None, **kwargs):
        """
        
        Initialises the MiriDroopModel class.
        
        Parameters: See class doc string.

        """
        super(MiriDroopModel, self).__init__(init=init, **kwargs)

        # Data type is droop.
        self.meta.model_type = 'DROOP'
        self.meta.reftype = 'DROOP'
        
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

        if droop_table is not None:
            try:
                self.droop_table = droop_table
            except (ValueError, TypeError) as e:
                strg = "droop_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
#         
#         # Copy the table column units, if defined.
#         droop_units = self.set_table_units('droop_table')
        
    def __str__(self):
        """
        
        Return the contents of the droop object as a readable
        string.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')

        # Describe the droop table
        if self.droop_table is not None:
            strg += self.get_data_str('droop_table', underline=True, underchar="-")
        else:
            strg += "No droop table!"
        return strg


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriDroopModel module.")
    
    PLOTTING = False
    SAVE_FILES = False

    droopdata = [(0.012,  0.012)]

    print("\nDroop calibration with factors derived from list of tuples:")
    with MiriDroopModel( droop_table=droopdata, detector='MIRIFUSHORT' ) as testdroop1:
        print(testdroop1)
        schema_names = testdroop1.get_field_names('droop_table')
        print("Field names in schema are:", schema_names)
        if PLOTTING:
            testdroop1.plot(description="testdroop1")
        if SAVE_FILES:
            testdroop1.save("test_droop_model1.fits", overwrite=True)
        del testdroop1
        
    print("Test finished.")
