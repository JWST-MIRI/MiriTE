#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An experimental extension to the standard STScI data model, which 
defines MIRI ramp jump correction models.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

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

@author: Vincent Geers (DIAS)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

#import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
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
        self.meta.model_type = 'JUMP'
        self.meta.reftype = 'JUMP'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'
        
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
        
    def __str__(self):
        """
        
        Display the contents of the jump object as a readable
        string.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')

        # Describe the fine jump table
        strg += "\nColumn names: " + str(self.fieldnames) + "\n"
        if self.finejump_table is not None:
            strg += self.get_data_str('finejump_table', underline=True, underchar="-")
        return strg


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
