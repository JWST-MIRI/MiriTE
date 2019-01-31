#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which 
defines MIRI latent decay models.

:Reference:

The STScI jwst.datamodels documentation.
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

27 Feb 2013: Created
01 Mar 2013: Reformatted the table and fixed parameters moved to metadata.
17 May 2013: Do not allow a blank table to be created.
22 Aug 2013: columnnames renamed to fieldnames.
02 Sep 2013: Pass the responsibility for creating record arrays to jwst.datamodels
             - a solution to the "Types in column 0 do not match" problem
             suggested by Michael Droettboom at STScI.
12 Sep 2013: Make the creation of an empty table a warning rather than an
             error.
24 Sep 2013: Reformatted model to include 4 latent decay tables instead of 1.
25 Sep 2013: Modified the constructor to accept a list of table objects.
             Corrected some typos and misunderstandings.
27 Sep 2013: Restored the LTUSED flag (use_flag) back to being Y/N.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: TYPE and REFTYPE are no longer identical.
09 Jul 2015: Removed duplication of table units between schema and metadata.
             Units are now only defined in the metadata.
             Use of the fieldnames class variable removed from the code and
             deprecated. It is now used only by a few conversion scripts.
11 Sep 2015: Removed duplicated plot method.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
14 Nov 2018: Explicitly set table column units based on the tunit definitions
             in the schema.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).

@author: Steven Beard (UKATC), Vincent Geers (DIAS)

"""
# This module is now converted to Python 3.


import warnings
#import numpy as np
#import numpy.ma as ma

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriLatentDecayModel']

class MiriLatentTable(object):
    """
    
    A helper class for generating a description of a MIRI latent table.
    The class may be used to construct a MiriLatentDecayModel object
    from a list of related tables.
    
    All the constructor parameters are optional. Any parameter not
    specified here will be given a None value, which will then be
    interpreted by the MiriLatentDecayModel class as a default value.
    
    :Parameters:
    
    regime: string, optional
        Description of the latent decay regime described by the table
        (e.g. "saturated", "unsaturated", etc...).
    exposure: int, optional
        ID of the exposure analyzed (e.g. 12029).
    background: float, optional
        Background level for latent decay measurement (in the units
        expected by the latent decay model - e.g. electrons/s/pixel)
    latent: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (parameter:object, factor:number,
        uncertainty:number), giving the latent decay parameters.
        Or: A numpy record array containing the same information as above.
        Or: None if this table is undefined.

    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """

    def __init__(self, regime=None, exposure=None, background=None,
                 latent=None):
        """
        
        Initialises the MiriLatentTable class.
        
        Parameters: See class doc string.

        """
        self.regime = regime
        self.exposure = exposure
        self.background = background
        self.latent = latent
        

class MiriLatentDecayModel(MiriDataModel):
    """
    
    A generic data model for a MIRI latent decay table, based on the STScI
    base model,vDataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the latent parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    latent_tables: list of MiriLatentTable objects (optional)
        A list of latent table definitions, containing metadata
        parameters plus a decay table definition for each latent
        table defined in the schema.
        If fewer latent tables are defined here than in the schema,
        the extra tables in the schema are left undefined (or defined
        from the initializer).
        If more latent tables are defined here than in the schema,
        the extra tables defined here will be ignored (and a warning
        issued).       
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_latent.schema.yaml"
    fieldnames = ('TYPE', 'T_OFFSET', 'DECAY_PEAK', 'DECAY_FREQ', 'TAU')
   
    def __init__(self, init=None, latent_tables=None, **kwargs):
        """
        
        Initialises the MiriLatentDecayModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLatentDecayModel, self).__init__(init=init, **kwargs)

        # Data type is latent decay.
        self.meta.reftype = 'LATENT'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type:
            self.meta.model_type = model_type        

        # This is a reference data model.
        self._reference_model()

        # Fill in the latent tables from the information provided.
        if latent_tables is not None:
            if isinstance(latent_tables, (tuple,list)):
                tableno = 0
                for latent_table in latent_tables:
                    # Assume each item in the latent_tables list matches
                    # one of the latent tables defined in the schema, in
                    # the same order as defined: latent1, laten2, etc...
                    tableno += 1
                    tablename = "latent%d" % tableno
                    tablemetadata = "meta.latent%d" % tableno
                    
                    # Skip a definition if it is not provided. It is assumed
                    # the table will be defined somewhere else.
                    if latent_table is not None:
                        # Each definition, if it exists, must be a
                        # MiriLatentTable object.
                        if isinstance(latent_table, MiriLatentTable):
                            # Check that the data structure has a matching
                            # attribute.
                            if hasattr(self, tablename):
                                # There is a matching attribute. Try to
                                # write a new table.
                                try:
                                    self[tablename] = latent_table.latent
                                except (ValueError, TypeError) as e:
                                    strg = "latent table %s" % tablename
                                    strg +=" must be a numpy record array"
                                    strg += " or list of records."
                                    strg += "\n   %s" % str(e)
                                    raise TypeError(strg)
                                # If the table has been written successfully,
                                # update the metadata.
                                metaregime = tablemetadata + ".regime"
                                metaexposure = tablemetadata + ".exposure"
                                metabackground = tablemetadata + ".background"
                                metaused = tablemetadata + ".use_flag"
                                try:
                                    self[metaregime] = latent_table.regime
                                    self[metaexposure] = latent_table.exposure
                                    self[metabackground] = latent_table.background
                                    self[metaused] = 'Y'
                                except (ValueError, TypeError) as e:
                                    strg = "Could not define metadata for "
                                    strg = "latent table %s" % tablename
                                    strg += "\n   %s" % str(e)
                                    raise TypeError(strg)
                            else:
                                # There isn't a matching attribute.
                                strg = "\n***Data defined for latent table %s " % \
                                    tablename
                                strg += "ignored. Not defined in schema."
                                warnings.warn(strg)
                        else:
                            strg = "latent table %s should " % tablename
                            strg += "either be None or a MiriLatentTable. "
                            strg += "%s given instead." % \
                                latent_table.__class__.__name__
                            raise TypeError(strg)                            
            else:
                strg = "latent_tables must either be a list of table "
                strg += "definitions or None. %s given instead." % \
                    latent_tables.__class__.__name__
                raise TypeError(strg)

        # Copy the table column units and at the same time
        # check for a data object containing entirely empty tables.
        tableno = 1
        tablename = "latent%d" % tableno
        allempty = True
        while hasattr(self, tablename):
            latent = getattr(self, tablename)
            if latent is not None and len(latent) > 1:
                allempty = False
                # Copy the table column units from the schema, if defined.
                latent_units = self.set_table_units(tablename)
            tableno += 1
            tablename = "latent%d" % tableno
        if allempty:
            strg = "\n***MiriLatentDecayModel object created with all empty tables."
            warnings.warn(strg)
      
    def __str__(self):
        """
        
        Return the contents of the latent decay object as a readable
        string.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title_and_metadata()

        # Describe the latent decay tables
        tableno = 1
        tablename = "latent%d" % tableno
        while hasattr(self, tablename):
            if getattr(self,tablename) is not None:
                strg += self.get_data_str(tablename, underline=True, underchar="-")
            tableno += 1
            tablename = "latent%d" % tableno
        return strg

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriLatentDecayModel module.")
    
    PLOTTING = False
    SAVE_FILES = False
    
    regime1 = "Unsaturated"
    exposure_id1 = 12029
    background1 = 8.34446
    latent_params1 = [('Fast',        355, 25000, 0.12,   8.33),
                     ('Intermediate', 365,    35, 0.01,   100),
                     ('Slow',         380,     9, 0.0029, 345)
                     ]
    latent_table1 = MiriLatentTable( regime1, exposure_id1, background1,
                                    latent_params1 )

    regime2 = "Full well"
    exposure_id2 = 12030
    background2 = 8.34446
    latent_params2 = [('Fast',        370, 28000, 0.12,   8.33),
                     ('Intermediate', 380,    16, 0.01,   100),
                     ('Slow',         380,    12, 0.0029, 345)
                     ]
    latent_table2 = MiriLatentTable( regime2, exposure_id2, background2,
                                    latent_params2 )
    latent_tables = [latent_table1, latent_table2] 

    print("\nLatent model 1: Two tables")
    with MiriLatentDecayModel( latent_tables=latent_tables) as testlatent1:
        print(testlatent1)
        if PLOTTING:
            testlatent1.plot(description="testlatent1")
        if SAVE_FILES:
            testlatent1.save("test_latent_model1.fits", overwrite=True)
        del testlatent1

    print("\nLatent model 2: No tables")
    with MiriLatentDecayModel( latent_tables=None) as testlatent2:
        print(testlatent2)
        if PLOTTING:
            testlatent2.plot(description="testlatent2")
        if SAVE_FILES:
            testlatent2.save("test_latent_model2.fits", overwrite=True)
        del testlatent2

    print("\nLatent model 3: Five tables")
    latent_tables = [latent_table1, latent_table2, latent_table1,
                     latent_table2, latent_table1] 
    with MiriLatentDecayModel( latent_tables=latent_tables) as testlatent3:
        print(testlatent3)
        if PLOTTING:
            testlatent3.plot(description="testlatent3")
        if SAVE_FILES:
            testlatent3.save("test_latent_model3.fits", overwrite=True)
        del testlatent3

    print("Test finished.")
