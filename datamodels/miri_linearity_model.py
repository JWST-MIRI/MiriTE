#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI linearity 
correction data.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

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
13 Feb 2018: Corrected a typo in the coeffs.setter. Added plot_goeffs,
             get_coeffs_str, apply, get_forward_table and get_reverse_table
             functions and updated the tests run in the main program.
20 Apr 2018: Corrected a bug where the second to last entry of the reverse_table
             could end up containing zero.
18 May 2018: Use explicit rounding when converting floating point to integer
             in get_forwards_table and get_reverse_table.
06 Jul 2018: Merged schema with JWST software.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).
04 Oct 2019: Added missing primary array name.
07 Oct 2019: Updated linearity_reference_flags to include only standard flag
             names. Removed '.yaml' suffix from schema references.

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""

import numpy as np
import numpy.ma as ma
import sys

# Import the MIRI measured data model.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_measured_model import MiriMeasuredModel

# List all classes and global functions here.
__all__ = ['linearity_reference_flags', 'MiriLinearityModel']

# The new JWST linearity reference flags
linearity_reference_setup = \
            [(0, 'DO_NOT_USE',      'Bad pixel. Do not use.'),
             (1, 'NONLINEAR',       'Pixel highly non-linear'),
             (2, 'NO_LIN_CORR',     'Linearity correction not available'),
             (3, 'OTHER_BAD_PIXEL', 'Bad solution'),                       # Was CDP_BAD_SOLUTION
             (4, 'DROPOUT',         'Data derived from incomplete input')] # Was CDP_PARTIAL_DATA 
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
    schema_url = "miri_linearity.schema"
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
        self.meta.reftype = 'LINEARITY'
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

        # This is a reference data model.
        self._reference_model()
                    
        if fitref is not None:
            self.meta.fit.reference = fitref
        if fitmodel is not None:
            self.meta.fit.model = fitmodel
        if order is not None:
            self.meta.fit.order = order

    def get_primary_array_name(self):
        """

        Returns the name "primary" array for this model.

        """
        return 'data'


    def __str__(self):
        """
        
        Return the contents of the linearity reference object as a readable
        string.
        
        :Returns:
        
        description: str
            A string containing a summary of the contents.
        
        """
        # First obtain a string describing the underlying measured
        # model.
        strg = super(MiriLinearityModel, self).__str__()
        
        # Add the extras
        if (self.data is not None) and (self.data.shape[2] > 2):
            strg += self.get_coeffs_str(underline=True)
            strg += "\n"
        if self.meta.fit.model is not None and self.meta.fit.order is not None:
            strg += "\nFit model is \'%s\' of order %d.\n" % \
                (self.meta.fit.model, self.meta.fit.order)
        elif self.meta.fit.model is not None:
            strg += "\nFit model is \'%s\' of UNKNOWN order.\n" % \
                self.meta.fit.model
        else:
            strg += "\nFit model is UNDEFINED.\n"
        return strg
    
    def plot_coeffs(self, row, column, description=""):
        """
        
        Plot a graph showing the shape of the linearity function at the
        specified row and column.
        
        :Parameters:
        
        row: int
            The row at which the polynomial is to be plotted.
        column: int or tuple of its.
            The column at which the polynomial is to be plotted.
        description: str (optional)
            An optional description to be appended to the plot title.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt
        
        # Start with an array of smoothly increasing DN values
        xarray = np.arange(0.0, 65000.0, 1000.0)
        
        # Convert the array using the given polynomial coefficients
        ncoeffs = self.data.shape[0]
        yarray = np.zeros_like(xarray)
        for icoeff in range(0, ncoeffs):
            yarray += self.data[icoeff,row,column] * (xarray ** icoeff)     

        # Define fixed labels
        xlabel = "Input DN"
        ylabel = 'Output DN'
        tstrg = \
            "Linearity polynomial for detector %s" % str(self.meta.instrument.detector) + \
            " (row=%d,col=%d)" % (row,column)
        if description:
            tstrg += "\n" + description
            
        # Plot the polynomial as an XY plot.
        mplt.plot_xy(xarray, yarray, linefmt='bo', xlabel=xlabel,
                     ylabel=ylabel, title=tstrg)
        mplt.show_plot()
        
    def get_coeffs_str(self, underline=False):
        """
        
        Return a string summarising the polynomial coefficients contained
        in the linearity data model.
        
        :Parameters:

        underline: bool, optional
            Set to True if the title should be underlined. The default
            is False.
        
        :Returns:
        
        coeffs_str: str
            A string summarising the polynomial coefficients contained.
        
        """
        ncoeffs = self.data.shape[0]        
        # Obtain the two arrays containing the median coefficients
        # for the left and right sides of the detector array
        mcoeffs_l = []
        mcoeffs_r = []
        boundary2 = self.data.shape[2]//2
        for icoeff in range(0, ncoeffs):
            mcoeffs_l.append(np.median(self.data[icoeff,:,:boundary2]))
            mcoeffs_r.append(np.median(self.data[icoeff,:,boundary2+1:]))
        
        # Determine whether the two arrays are identical or different
        identical = True
        for icoeff in range(0, ncoeffs):
            if abs(mcoeffs_l[icoeff] - mcoeffs_r[icoeff]) > sys.float_info.epsilon:
                identical = False
        
        strg = "\nMedian polynomial coefficients"
        if underline:
            strg += "\n-----------------------------"
            strg += "\n      "
        for icoeff in range(0, ncoeffs):
            strg += "         L%d  " % icoeff
        if identical:
            strg += "\nAll:  "
            for coeff in mcoeffs_l:
                strg += " %12.6g" % coeff
        else:
            strg += "\nLeft: "
            for coeff in mcoeffs_l:
                strg += " %12.6g" % coeff
            strg += "\nRight:"
            for coeff in mcoeffs_r:
                strg += " %12.6g" % coeff
        return strg
    
    def apply(self, inramp):
        """
        
        Apply the linearity coefficients to a ramp of DN values.
        The input data must match the number of rows and columns
        contained in the linearity coefficients array.
        
        :Parameters:
        
        inramp: array-like
            Input array containing an uncorrected ramp of DN values.
        
        :Returns:
        
        outramp: array-like
            Output array containing a corrected ramp of DN values.
        
        """
        inramp = np.asarray(inramp)
        assert (len(inramp.shape) > 1)
        assert (inramp.shape[-1] == self.data.shape[-1])
        assert (inramp.shape[-2] == self.data.shape[-2])
        ncoeffs = self.data.shape[0]
        outramp = np.zeros_like(inramp)
        for icoeff in range(0, ncoeffs):
            outramp += self.data[icoeff,:,:] * (inramp[:,:,:] ** icoeff)     
        return outramp
    
    def get_forward_table(self, row, column, max_dn=65535):
        """
        
        Return a table which translates input DN to output DN using
        the coefficients associated with the given row and column.
        The table is contained in an array of integers so that the
        expression
        
           output_dn = forward_table[ input_dn ]
           
        will correct the input DN value for nonlinearity.
        The forward table is useful in situations where it is quicker
        to calculate all possible polynomial conversions in advance
        (e.g. if the conversion is to be applied more than max_dn times).
        
        :Parameters:
        
        row: int
            The row at which the polynomial is to be plotted.
        column: int or tuple of its.
            The column at which the polynomial is to be plotted.
        max_dn: int (optional)
            The maximum input DN value. Also the size of the returned table.
            Defaults to 65525. The minimum DN value is always 0.
            
        :Returns:
        
        forward_table: array of int
            An array of integers containing the translation table.
            
        """        
        # Start with an array of all possible DN values
        inarray = np.arange(0.0, float(max_dn), 1.0)
        
        # Convert the array using the given polynomial coefficients
        ncoeffs = self.data.shape[0]
        farray = np.zeros_like(inarray)
        for icoeff in range(0, ncoeffs):
            farray += self.data[icoeff,row,column] * (inarray ** icoeff)
            
        # Convert the output array to integer.
        outarray = np.floor(farray+0.5).astype(np.int)
        return outarray

    def get_reverse_table(self, row, column, max_dn=65535, fill_gaps=True):
        """
        
        Return a table which translates output DN back to input DN using
        the coefficients associated with the given row and column.
        The table is contained in an array of integers so that the
        expression
        
           input_dn = forward_table[ output_dn ]
           
        will reveal the input DN value that would generate the given
        output when converted with the nonlinearity polynomial.
        The reverse table is useful in situations (such as in a data
        simulator) where it is necessary to run the pipeline steps
        backwards.
        
        :Parameters:
        
        row: int
            The row at which the polynomial is to be plotted.
        column: int or tuple of its.
            The column at which the polynomial is to be plotted.
        max_dn: int (optional)
            The maximum input DN value.
            Defaults to 65525. The minimum DN value is always 0.
        fill_gaps: bool (optional)
            If True (the default), interpolate the table to fill any gaps.
            
        :Returns:
        
        reverse_table: array of int
            An array of integers containing the translation table.
            
        """
        # Determine the maximum DN in the reverse table
        ncoeffs = self.data.shape[0]
        fmax = 0.0
        for icoeff in range(0, ncoeffs):
            fmax += self.data[icoeff,row,column] * (float(max_dn) ** icoeff)
        max_dn_out = int(fmax+0.5)
        # Initialise a reverse table
        rarray = np.zeros([max_dn_out+1])
        
        # Start by generating a forward table.
        ftable = self.get_forward_table(row, column, max_dn=max_dn)
        
        # Use the forward table to populate the reverse table
        for indn in range(0, max_dn):
            outdn = int(ftable[indn])
            rarray[outdn] = indn
            
        # If any elements remain unpopulated, interpolate the nearest
        # non-zero values.
        if fill_gaps:
            runopen = False
            for outdn in range(1, max_dn_out):
                if rarray[outdn] == 0:
                    if not runopen:
                        runstart = outdn
                        runopen = True
                else:
                    if runopen:
                        runend = outdn
                        runopen = False
                        leftvalue = rarray[runstart-1]
                        rightvalue = rarray[runend]
                        inc = (rightvalue - leftvalue)/float(runend-runstart)
                        for rundn in range(runstart,runend):
                            rarray[rundn] = rarray[rundn-1] + inc
            if rarray[max_dn_out] == 0:
                rarray[max_dn_out] = max_dn
            if rarray[max_dn_out-1] == 0:
                rarray[max_dn_out-1] = (rarray[max_dn_out-2] + rarray[max_dn_out])/2.0
            
        # Convert the output array to integer.
        reverse_array = np.floor(rarray+0.5).astype(np.int)
        return reverse_array

    # "coeffs" is an alias for the "data" attribute.
    @property
    def coeffs(self):
        if hasattr(self, 'data'):
            return self.data
        else:
            return None

    @coeffs.setter
    def coeffs(self, data):
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

    coeffs1d = [0.0, 0.865384, 4.64239e-6, -6.16093e-11, 7.23130e-16]
    
    coeffs3d = np.asarray([[coeffs1d, coeffs1d, coeffs1d],
                          [coeffs1d, coeffs1d, coeffs1d],
                          [coeffs1d, coeffs1d, coeffs1d],
                          ])
    coeffs3x3x5 = np.moveaxis(coeffs3d,-1,0)
    err3x3x5 = 0.01 * np.ones_like(coeffs3x3x5)
    dq3x3 = np.asarray([[0,1,0],[1,0,1],[0,1,0]])
    
    testdata = np.asarray([[19000.0, 20000.0, 21000.0],
                           [22000.0, 18000.0, 20000.0],
                           [22000.0, 18000.0, 20000.0]])
    testramp = np.asarray([0.2*testdata,
                           0.4*testdata, 0.6*testdata,
                           0.8*testdata,
                           testdata])

    print("Linearity data with data + err + dq:")
    with MiriLinearityModel( coeffs=coeffs3x3x5, err=err3x3x5, dq=dq3x3,
                                dq_def=linearity_reference_flags ) \
            as testdata1:
        print(testdata1)
        # Test the apply function
        print("input ramp=", testramp)
        outramp = testdata1.apply( testramp )
        print("output ramp=", outramp)
        # Test the calculation of forward and reverse tables
        forward = testdata1.get_forward_table(row=1, column=1)
        print("Forward table of length", len(forward))
        print(str(forward))
        reverse = testdata1.get_reverse_table(row=1, column=1)
        print("Reverse table of length", len(reverse))
        print(str(reverse))
        
        #print(str(reverse))
        if PLOTTING:
            testdata1.plot_coeffs(row=1, column=1, description="testdata1")
            testdata1.plot(description="testdata1")
        if SAVE_FILES:
            testdata1.save("test_linearity_model1.fits", overwrite=True)
        del testdata1

    if READ_CDP_DATA:
        from miri.datamodels.cdplib import get_cdp
        print("Reading linearity data from a CDP files")
        for detector in ('MIRIMAGE', 'MIRIFUSHORT', 'MIRIFULONG'):
            with get_cdp('LINEARITY', detector=detector) as testdata3:
                print(testdata3)
                # Test the calculation of forward and reverse tables
                forward = testdata3.get_forward_table(row=512, column=400)
                print("Forward table of length", len(forward))
                print(str(forward))
                reverse = testdata3.get_reverse_table(row=512, column=400)
                print("Reverse table of length", len(reverse))
                print(str(reverse))
                # Test the calculation of forward and reverse tables
                forward = testdata3.get_forward_table(row=512, column=650)
                print("Forward table of length", len(forward))
                print(str(forward))
                reverse = testdata3.get_reverse_table(row=512, column=650)
                print("Reverse table of length", len(reverse))
                print(str(reverse))
                if PLOTTING:
                    testdata3.plot_coeffs(row=512, column=400, description="testdata3")
                    testdata3.plot_coeffs(row=512, column=650, description="testdata3")
                    testdata3.plot(description="testdata3")
    
#                 print("All occurences of BUNIT: " + str(testdata3.find_fits_values('BUNIT')))
#                 print("Units of the COEFFS data: " + str(testdata3.get_fits_keyword('BUNIT', 'COEFFS')))
#     
#                 print("All COMMENT records: " + testdata3.get_comments_str())
#                 print("All HISTORY records: " + testdata3.get_history_str())
#                 print("Getting ORDER keyword: " + str(testdata3.get_fits_keyword('ORDER')))
#                 testdata3.set_fits_keyword('ORDER', 3)
#                 print("Getting ORDER keyword again: " + str(testdata3.get_fits_keyword('ORDER')))
#                 testdata3.add_history('Mucked about with during testing.')
#                 print("All HISTORY records again: " + testdata3.get_history_str())
            
                del testdata3

    print("Test finished.")
