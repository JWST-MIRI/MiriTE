#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI linearity 
correction data.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

10 Jan 2013: Created
21 Jan 2013: Included data type and fit parameters in metadata.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
26 Feb 2013: BIAS image removed on advice from Jane Morrison.
             Default fit model and order should be None.
01 Mar 2013: Updated to match Jane Morrison's latest DHAS flags.
07 Jun 2013: Extra optional tests with CDP data added at the end.
04 Jul 2013: Don't display masked data when it is the same as
             unmasked data. 
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
10 Apr 2014: Modified for jsonschema draft 4: Functions made more
             independent of schema structure. Modified to define
             data units using the set_data_units method.
09 Jul 2014: JWST reference flags table added.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out. Corrected REFTYPE.
10 Oct 2014: Expanded the description of FITMODEL and added FITREF for
             consistency with the distortion models.
29 Oct 2014: Disabled FITERR, MIN and MAX from being automatically included
             in model. Renamed SCI plane to COEFFs.
04 Nov 2014: Added coeffs as an alias for data.
09 Dec 2014: Obsolete field_def table removed.
10 Jun 2015: Obsolete calls to fits_history and fits_comment functions
             replaced by calls to generic functions.
20 Aug 2015: Name changed from MiriNonlinearityModel to MiriLinearityModel,
             for consistency with STScI naming. minimage, maximage and fiterr
             arrays removed.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
20 Nov 2017: Updated the schema so the world coordinates are written to
             the COEFFS data component.

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import numpy as np
import numpy.ma as ma

# Import the MIRI measured data model.
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_measured_model import MiriMeasuredModel

# List all classes and global functions here.
__all__ = ['linearity_reference_flags', 'MiriLinearityModel']

# The new JWST linearity reference flags
linearity_reference_setup = \
            [(0, 'DO_NOT_USE',  'Bad pixel. Do not use.'),
             (1, 'NONLINEAR',   'Pixel highly non-linear'),
             (2, 'NO_LIN_CORR', 'Linearity correction not available'),
             (3, 'CDP_BAD_SOLUTION', 'Bad solution'),
             (4, 'CDP_PARTIAL_DATA', 'Data derived from incomplete input')]
linearity_reference_flags = insert_value_column( linearity_reference_setup )


class MiriLinearityModel(MiriMeasuredModel):
    """
    
    A data model for MIRI linearity correction data.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the coeffs parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    coeffs: numpy array (optional)
        A 3-D array containing the linearity coefficients.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    err: numpy array (optional)
        An array containing the uncertainty in the the linearity
        coefficients.
        Must be broadcastable onto the coeffs array.
    dq: numpy array (optional)
        An array containing the quality of the linearity coefficients.
        Must be broadcastable onto the coeffs array.
    dq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, title:str),
        giving the meaning of values stored in the data quality array. For
        example: [(0, 'good','Good data'), (1, 'dead', 'Dead Pixel'),
        (2, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags.
    fitref: str (optional)
        A string containing a human-readable reference to a document
        describing the distortion model.
    fitmodel: str (optional)
        If a recognised JWST fitting model has been used (e.g. one of the
        models in the astropy.modeling package) a unique, machine-readable
        string defining the model used. If the model is not known or
        doesn't match one of the standard JWST models, leave this keyword
        blank and describe the model using the fitref parameter (above).
        Some example strings from astropy.modeling: Chebyshev1D', 'Chebyshev2D',
        'InverseSIP', 'Legendre1D','Legendre2D', 'Polynomial1D',
        'Polynomial2D', etc...
    order: int: optional
        The order of fit.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_linearity.schema.yaml"
    _default_dq_def = linearity_reference_flags

    def __init__(self, init=None, coeffs=None, dq=None, err=None,
                 dq_def=None, fitref=None, fitmodel=None,
                 order=None, **kwargs):
        """
        
        Initialises the MiriLinearityModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLinearityModel, self).__init__(init=init, data=coeffs,
                                                 dq=dq, err=err, dq_def=dq_def,
                                                 **kwargs)
        
        # Data type is non-linearity.
        self.meta.model_type = 'LINEARITY'
        self.meta.reftype = 'LINEARITY'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'
                    
        if fitref is not None:
            self.meta.fit.reference = fitref
        if fitmodel is not None:
            self.meta.fit.model = fitmodel
        if order is not None:
            self.meta.fit.order = order

    def __str__(self):
        """
        
        Display the contents of the dark reference object as a readable
        string.
        
        """
        # First obtain a string describing the underlying measured
        # model.
        strg = super(MiriLinearityModel, self).__str__()
        
        # Add the extras
        if self.meta.fit.model is not None and self.meta.fit.order is not None:
            strg += "\nFit model is \'%s\' of order %d.\n" % \
                (self.meta.fit.model, self.meta.fit.order)
        elif self.meta.fit.model is not None:
            strg += "\nFit model is \'%s\' of UNKNOWN order.\n" % \
                self.meta.fit.model
        else:
            strg += "\nFit model is UNDEFINED.\n"
        return strg

    # "coeffs" is an alias for the "data" attribute.
    @property
    def coeffs(self):
        if hasattr(self, 'data'):
            return self.data
        else:
            return None

    @coeffs.setter
    def mask(self, data):
        self.data = data


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriLinearityModel module.")

    PLOTTING = False
    SAVE_FILES = False
    READ_CDP_DATA = False

    data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])
    data3x3x2 = [data3x3,data3x3]

    print("Linearity data with data + err + dq:")
    with MiriLinearityModel( coeffs=data3x3x2, err=err3x3, dq=dq3x3,
                                dq_def=linearity_reference_flags ) \
            as testdata1:
        print(testdata1)
        if PLOTTING:
            testdata1.plot(description="testdata1")
        if SAVE_FILES:
            testdata1.save("test_linearity_model1.fits", overwrite=True)
        del testdata1

    if READ_CDP_DATA:
        print("Reading linearity data from a CDP file")
        with MiriLinearityModel( init="D:/MIRI-CDP-4/MIRI_FM_MIRIFUSHORT_12_LINEARITY_04.02.00.fits") \
                as testdata3:
            print(testdata3.get_meta_str())
            print("All occurences of BUNIT: " + str(testdata3.find_fits_values('BUNIT')))
            print("Units of the COEFFS data: " + str(testdata3.get_fits_keyword('BUNIT', 'COEFFS')))

            print("All COMMENT records: " + testdata3.get_comments_str())
            print("All HISTORY records: " + testdata3.get_history_str())
            print("Getting ORDER keyword: " + str(testdata3.get_fits_keyword('ORDER')))
            testdata3.set_fits_keyword('ORDER', 3)
            print("Getting ORDER keyword again: " + str(testdata3.get_fits_keyword('ORDER')))
            print("Getting ORDER keyword again: " + str(testdata3.get_fits_keyword('ORDER')))
        
            testdata3.add_fits_comment('Test comment added during development')
            testdata3.add_fits_comment('Another comment added during development')
            testdata3.add_fits_comment('The science data are cool.', hdu_name='COEFFS')
            testdata3.add_fits_comment('There is no data quality data.', hdu_name='DQ')
            print("All COMMENT records again: " + testdata3.get_comments_str())

            testdata3.add_history('Mucked about with during testing.')
            print("All HISTORY records again: " + testdata3.get_history_str())
        
            #print(str(testdata3.to_flat_dict()))
            #testdata3.search_schema('order', verbose=True)
        
            del testdata3

    print("Test finished.")
