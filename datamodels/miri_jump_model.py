#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An experimental extension to the standard STScI data model, which 
defines MIRI ramp jump correction models.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

25 Sep 2013: Created
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: TYPE and REFTYPE are no longer identical.
09 Dec 2014: Removed old commented out warnings.
11 Sep 2015: Removed duplicated plot method.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
14 Nov 2018: Explicitly set table column units based on the tunit definitions
             in the schema. Removed redundant function.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).

@author: Vincent Geers (DIAS)

"""

#import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriJumpModel']


class MiriJumpModel(MiriDataModel):
    """
    
    A generic data model for a MIRI ramp jump threshold table, based on 
    the STScI base model, DataModel.
    
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
        
    finejump_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (method:string, 
        rej_thresh:number), giving the rejection thresholds valid for different
        methods of identifying ramp jumps.
        Or: A numpy record array containing the same information as above.
        A finejump table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
    crthresh: number (optional)
        The count rate threshold in [e-/s] at which to switch from one method
        to the other.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_jump.schema.yaml"
    fieldnames = ('METHOD', 'REJ_THRESH')
    
    def __init__(self, init=None, finejump_table=None, crthresh=None, **kwargs):
        """
        
        Initialises the MiriJumpModel class.
        
        Parameters: See class doc string.

        """
        super(MiriJumpModel, self).__init__(init=init, **kwargs)

        # Data type is jump.
        self.meta.reftype = 'JUMP'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type is not None:
            self.meta.model_type = model_type        

        # This is a reference data model.
        self._reference_model()
        
        # Define the crthresh identifier, if specified. N.B. this ID is 
        # compulsory, so it must be specified either here, in the
        # source file or later after creation of this data object.
        if crthresh is not None:
            self.meta.finejump_table.crthresh = crthresh

        if finejump_table is not None:
            try:
                self.finejump_table = finejump_table
            except (ValueError, TypeError) as e:
                strg = "finejump_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
            
        # Copy the table column units from the schema, if defined.
        finejump_table_units = self.set_table_units('finejump_table')


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriJumpModel module.")
    
    PLOTTING = False
    SAVE_FILES = False

    finejumpdata = [('2PNTDIFF', 3.5),
                    ('YINT', 3.5)]
    crthresh = 5.0

    with MiriJumpModel( finejump_table=finejumpdata, crthresh = crthresh ) as testjump1:
        print(testjump1)
        if PLOTTING:
            testjump1.plot(description="testjump1")
        if SAVE_FILES:
            testjump1.save("test_jump_model1.fits", overwrite=True)
        del testjump1
        
    print("Test finished.")
