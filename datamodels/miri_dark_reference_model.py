#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI dark reference 
data. Essentially the same as the MIRI measured image data model, with 
additional metadata.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

19 Dec 2012: Created
21 Jan 2013: Included data type in metadata.
23 Jan 2013: Added plotting.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
21 Feb 2013: Base a dark reference model on MiriMeasuredModel rather than
             MiriRampModel, to allow for 3-D dark data.
04 Mar 2013: Only check the dimensionality of the main data array if it
             has been explicitly provided.
04 Jul 2013: Don't display masked data when it is the same as
             unmasked data. 
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
10 Apr 2014: Modified for jsonschema draft 4: Functions made more
             independent of schema structure. Modified to define data
             units using the set_data_units method.
09 Jul 2014: JWST reference flags table added.
29 Aug 2014: JWST reference flags table updated.
             Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out.
28 Oct 2014: Disabled automatic inclusion of the FITERR extension.
16 Jan 2015: Ensure the main data array is 3-D or 4-D when specified
             explicitly or when read from a file.
20 Aug 2015: Duplicated parts of schema now reference STScI model.
11 Sep 2015: Removed duplicated plot method.
01 Jun 2016: DARK data can now be 2-D, 3-D or 4-D.
20 Jun 2016: DARK DQ array can also now be 2-D, 3-D or 4-D.
15 Jul 2016: Removed obsolete flags.
17 Mar 2017: Corrected some documentation typos.
15 Jun 2017: Added average_group and extract_group methods.
             Added a comment to the data model created by the
             average_group and extract_group methods.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
             Added average_groups function.
07 Jul 2017: Added averaged flag.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
06 Jul 2018: Merged schema with JWST software. DARK data is now only
             accepted with 4-D data, err and dq arrays.

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""
# This module is now converted to Python 3.


import numpy as np
#import numpy.ma as ma
#import numpy.linalg as LA

# Import the MIRI ramp data model utilities.
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_measured_model import MiriMeasuredModel

# List all classes and global functions here.
__all__ = ['dark_reference_flags', 'MiriDarkReferenceModel']

# The new JWST dark reference flags. Note that these flags are a subset
# of the JWST master flags table. The names and descriptions are the same
# but the bit values can be different.
dark_reference_setup = \
            [(0, 'DO_NOT_USE',       'Bad pixel. Do not use.'),
             (1, 'UNRELIABLE_DARK',  'Dark variance large')]
dark_reference_flags = insert_value_column(dark_reference_setup)

#
# Global helper function
#
def linear_regression(x, y):
    """
    
    Linear regression function. Fit an intercept and a slope
    to a set of x,y coordinates
    
    """
    #print("Fitting", y, "\nagainst", x)
    matrix = np.vstack( [x, np.ones_like(x)] ).T
    slope, intercept = np.linalg.lstsq(matrix,y)[0]
    #print("gives slope=", slope, "intercept=", intercept)
    return (slope, intercept)


class MiriDarkReferenceModel(MiriMeasuredModel):
    """
    
    A data model for MIRI dark reference data.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the data parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    data: numpy array (optional)
        A 3-D or 4-D array containing the dark data.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    err: numpy array (optional)
        An array containing the uncertainty in the dark data.
        Must be broadcastable onto the data array.
    dq: numpy array (optional)
        An array containing the quality of the dark data.
        Must be broadcastable onto the data array.
    dq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, title:str),
        giving the meaning of values stored in the data quality array. For
        example: [(0, 'good','Good data'), (1, 'dead', 'Dead Pixel'),
        (2, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags.
    integration: number, optional
        The integration number for which this dark is valid.
    fitted_after_frame: number, optional
        The frame number beyond which the dark is fitted.
    averaged: bool, optional
        Indicates if the data model contains averaged dark data.
        The default is False.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_dark_reference.schema.yaml"
    _default_dq_def = dark_reference_flags

    def __init__(self, init=None, data=None, dq=None, err=None,
                 dq_def=None, integration=None, fitted_after=None,
                 averaged=False, **kwargs):
        """
        
        Initialises the MiriDarkReferenceModel class.
        
        Parameters: See class doc string.

        """
        super(MiriDarkReferenceModel, self).__init__(init=init, data=data,
                                                    dq=dq, err=err,
                                                    dq_def=dq_def,
                                                    **kwargs)

        # Data type is dark.
        if averaged:
            self.meta.model_type = 'DARK (averaged)'
        else:
            self.meta.model_type = 'DARK'
        self.meta.reftype = 'DARK'
        self.averaged = averaged
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'

        # Define the exposure type (if not already contained in the data model)
        if not self.meta.exposure.type:
            self.set_exposure_type( datatype='DARK' )
        
        if integration is not None:
            self.meta.integration_number = integration
        if fitted_after is not None:
            self.meta.fitted_after_frame = fitted_after

        # The main data array should be 2-D, 3-D or 4-D.
        # TODO: Can this check be included in the schema?
        if data is not None:
            if not hasattr(data, 'ndim'):
                data = np.asarray(data)
            if data.ndim < 2 or data.ndim > 4:
                strg = "The main data array in a dark reference object must be "
                strg += "2-D, 3-D or 4-D. %d-D data provided" % data.ndim
                raise TypeError(strg)
        elif self.data is not None and len(self.data) > 0 and \
             hasattr(self.data, 'ndim'):
            if self.data.ndim < 2 or self.data.ndim > 4:
                strg = "The main data array in a dark reference object must be "
                strg += "2-D, 3-D or 4-D. %d-D data provided" % self.data.ndim
                raise TypeError(strg)

#        if fiterr is not None:
#            self.fiterr = fiterr
#        self._fiterr_mask = None
#        self._fiterr_fill = 'max'
#        self._fiterr_fill_value = None
#
#        # Set the units of the fiterr array, if defined in the schema.
#        fitunits = self.set_data_units('fiterr')
        
    def __str__(self):
        """
        
        Return the contents of the dark reference object as a readable
        string.
        
        """
        # First obtain a string describing the underlying measured
        # model.
        strg = super(MiriDarkReferenceModel, self).__str__()
        
        # Add the extras
#        if self.fiterr is not None:
#            strg += self.get_data_str('fiterr', underline=True, underchar="-")
#            if self.maskable():
#                fiterr_masked = self.fiterr_masked
#                if fiterr_masked is not None:
#                    title = self.get_data_title('fiterr') + " (masked)"
#                    len2 = len(title)
#                    title += "\n" + (len2 * '~')
#                    strg += title + "\n" + str(fiterr_masked) + "\n"
        return strg

    def slope_data(self, startgroup=None, endgroup=None):
        """
        
        Return a new data model in which the group planes have
        been fitted with a straight line to make slope data.
        The units of the output data change from 'DN' to 'DN/s'

        NOTE: This function is very inefficient!
        
        :Parameters:
        
        startgroup: int, optional
            The first group to be averaged. Defaults to the first group
            in the data model.
        endgroup: int, optional
            The last group to be averaged. Defaults to the last group in
            the data model.
            
        :Returned:
        
        slope_model: MiriDarkReferenceModel
            A MiriDarkReferenceModel with a slope fitted to all
            the group data.
            
        """
        # First determine the group time from the detector readout mode
        grptime = 1.0
        if hasattr(self, 'meta') and hasattr( self.meta, 'exposure'):
            if hasattr(self.meta.exposure, 'group_time') and \
               self.meta.exposure.group_time is not None:
                grptime = self.meta.exposure.group_time
            elif hasattr( self.meta.exposure, 'readpatt') and \
                 self.meta.exposure.readpatt is not None:
                readpatt = self.meta.exposure.readpatt
                # TODO: Subarray modes not implemented yet
                if 'FAST' in readpatt:
                    grptime = 2.785
                else:
                    grptime = 27.85
            if hasattr(self.meta.exposure, 'group_time'):
                self.meta.exposure.group_time = grptime
        
        if self.data.ndim == 4:
            nints = self.data.shape[0]
            ngroups = self.data.shape[1]
            nrows = self.data.shape[2]
            ncolumns = self.data.shape[3]

            # Initialise the output.
            output = np.zeros([nints, nrows, ncolumns], dtype=self.data.dtype)

            # Full straight line fit
            if startgroup is None and endgroup is None:
                timearray = grptime * np.array( list(range(0, ngroups)) )
                for intg in range(0, nints):
                    for row in range(0, nrows):
                        for column in range(0, ncolumns):
                            # Ouch! Slope calculated one (row,column) at a time.
                            # Can the efficiency be improved?
                            (slope, ic) = linear_regression( timearray,
                                                self.data[intg,:,row,column] )
                            output[intg,row,column] = slope
            else:
                if startgroup is None:
                    startgroup = 0
                if endgroup is None:
                    endgroup = self.data.shape[1]
                timearray = grptime * np.array( list(range(startgroup, endgroup)) )
                for intg in range(0, nints):
                    for row in range(0, nrows):
                        for column in range(0, ncolumns):
                            # Ouch! Slope calculated one (row,column) at a time.
                            # Can the efficiency be improved?
                            (slope, ic) = linear_regression( timearray,
                                                self.data[intg,startgroup:endgroup,row,column] )
                            output[intg,row,column] = slope
            output_dq = np.squeeze(self.dq)
        
            # Create a new 3-D dark data model from this slope data.
            newmodel = MiriDarkReferenceModel( data=output, err=None,
                                               dq=output_dq )
            # Copy the metadata to the new data model.
            # NOTE: The WCS metadata may be copied incorrectly.
            newmodel.copy_metadata( self )
            if grptime == 1.0:
                newmodel.meta.data.units = 'DN/group'
            else:
                newmodel.meta.data.units = 'DN/s'
            newmodel.add_comment( "TGROUP assumed %.3fs" % grptime )
            historystrg = "New DARK made by calculating the slope for %d groups" % \
                ngroups
            newmodel.add_comment( historystrg )
            return newmodel
        
        elif self.data.ndim == 3:
            ngroups = self.data.shape[0]
            nrows = self.data.shape[1]
            ncolumns = self.data.shape[2]

            # Initialise the output.
            output = np.zeros([nrows, ncolumns], dtype=self.data.dtype)

            # Full straight line fit
            if startgroup is None and endgroup is None:
                timearray = grptime * np.array( list(range(0, ngroups)) )
                for row in range(0, nrows):
                    for column in range(0, ncolumns):
                        # Ouch! Slope calculated one (row,column) at a time.
                        # Can the efficiency be improved?
                        (slope, ic) = linear_regression( timearray,
                                            self.data[intg,:,row,column] )
                        output[intg,row,column] = slope
            else:
                if startgroup is None:
                    startgroup = 0
                if endgroup is None:
                    endgroup = self.data.shape[1]
                timearray = grptime * np.array( list(range(startgroup, endgroup)) )
                for row in range(0, nrows):
                    for column in range(0, ncolumns):
                        # Ouch! Slope calculated one (row,column) at a time.
                        # Can the efficiency be improved?
                        (slope, ic) = linear_regression( timearray,
                                            self.data[intg,startgroup:endgroup,row,column] )
                        output[intg,row,column] = slope
            output_dq = np.squeeze(self.dq)
        
            # Create a new 3-D dark data model from this slope data.
            newmodel = MiriDarkReferenceModel( data=output, err=None,
                                               dq=output_dq, averaged=True )
            # Copy the metadata to the new data model.
            # NOTE: The WCS metadata may be copied incorrectly.
            newmodel.copy_metadata( self )
            if grptime == 1.0:
                newmodel.meta.data.units = 'DN/group'
            else:
                newmodel.meta.data.units = 'DN/s'
            newmodel.add_comment( "TGROUP assumed %.3fs" % grptime )
            historystrg = "New DARK made by calculating the slope for %d groups" % \
                ngroups
            newmodel.add_comment( historystrg )
            return newmodel
        
        else:
            raise TypeError("Data model has wrong number of dimensions")

    def average_groups(self, startgroup=None, endgroup=None, normalize=False):
        """
        
        Return a new data model in which the groups have been averaged.
        A 4-D data model [integrations, groups, rows, columns] will be
        reduced to 3-D [integrations, rows, columns] and a 3-D model
        [groups, rows, columns] will be reduced to 2-D [rows, columns].
        
        :Parameters:
        
        startgroup: int, optional
            The first group to be averaged. Defaults to the first group
            in the data model.
        endgroup: int, optional
            The last group to be averaged. Defaults to the last group in
            the data model.
        normalize: bool, optional
            If True, normalize the data so the average is 1.0.
            The default is False.
            
        :Returned:
        
        averaged_model: MiriDarkReferenceModel
            A MiriDarkReferenceModel with all the group data averaged together.
            
        """
        if self.data.ndim == 4:
        
            # Determine a new 3-D dark data model by averaging together
            # the specified groups.
            if startgroup is None and endgroup is None:
                ngroups = self.data.shape[1]
                cube_data = np.sum(self.data, 1)
                cube_data = cube_data / float(ngroups)
                #cube_err = LA.norm(self.err, ord=2, axis=1)
                errsq = self.err * self.err
                cube_err = np.sum(errsq, 1)
                cube_err = np.sqrt( cube_err ) / float(ngroups)
            else:
                if startgroup is None:
                    startgroup = 0
                if endgroup is None:
                    endgroup = self.data.shape[1]
                ngroups = endgroup - startgroup + 1
                cube_data = np.sum(self.data[:,startgroup:endgroup,:,:], 1)
                cube_data = cube_data / float(ngroups)
                #cube_err = LA.norm(self.err, ord=2, axis=1)
                errsq = self.err[:,startgroup:endgroup,:,:] * \
                        self.err[:,startgroup:endgroup,:,:]
                cube_err = np.sum(errsq, 1)
                cube_err = np.sqrt( cube_err ) / float(ngroups)
            if normalize:
                posdata = np.where(cube_data > 0.0)
                factor = np.mean(cube_data[posdata])
#                 print("Normalizing by a factor of %f" % factor)
                if factor > 0.0:
                    cube_data = cube_data / factor
                    cube_err = cube_err / factor
            cube_dq = np.squeeze(self.dq)
            newmodel = MiriDarkReferenceModel( data=cube_data, err=cube_err,
                                               dq=cube_dq )
            # Copy the metadata to the new data model.
            # NOTE: The WCS metadata may be copied incorrectly.
            newmodel.copy_metadata( self )
            historystrg = "New DARK made by averaging groups %d to %d" % \
                (startgroup, endgroup)
            if normalize:
                historystrg += " (normalized to 1.0)"
            newmodel.add_comment( historystrg )
            return newmodel
        
        elif self.data.ndim == 3:
        
            # Determine a new 2-D dark data model by averaging together
            # all the groups.
            if startgroup is None and endgroup is None:
                ngroups = self.data.shape[0]
                image_data = np.sum(self.data, 0)
                image_data = image_data / float(ngroups)
                errsq = self.err * self.err
                image_err = np.sum(errsq, 0)
                image_err = np.sqrt( image_err ) / float(ngroups)
            else:
                if startgroup is None:
                    startgroup = 0
                if endgroup is None:
                    endgroup = self.data.shape[0]
                ngroups = endgroup - startgroup + 1
                image_data = np.sum(self.data[startgroup:endgroup,:,:], 0)
                image_data = image_data / float(ngroups)
                errsq = self.err[startgroup:endgroup,:,:] * \
                        self.err[startgroup:endgroup,:,:]
                image_err = np.sum(errsq, 0)
                image_err = np.sqrt( image_err ) / float(ngroups)
            if normalize:
                posdata = np.where(image_data > 0.0)
                factor = np.mean(image_data[posdata])
#                 print("Normalizing by a factor of %f" % factor)
                if factor > 0.0:
                    image_data = image_data / factor
                    image_err = image_err / factor
            image_dq = np.squeeze(self.dq)
            newmodel = MiriDarkReferenceModel( data=image_data, err=image_err,
                                               dq=image_dq, averaged=True )
            # Copy the metadata to the new data model.
            # NOTE: The WCS metadata may be copied incorrectly.
            newmodel.copy_metadata( self )
            historystrg = "New DARK made by averaging groups %d to %d" % \
                (startgroup, endgroup)
            if normalize:
                historystrg += " (normalized to 1.0)"
            newmodel.add_comment( historystrg )
            return newmodel
        
        else:
            raise TypeError("Data model has wrong number of dimensions")

    def extract_group(self, group=0, normalize=False):
        """
        
        Return a new data model containing a single group extracted
        from the current data model.
        A 4-D data model [integrations, groups, rows, columns] will be
        reduced to 3-D [integrations, rows, columns] and a 3-D model
        [groups, rows, columns] will be reduced to 2-D [rows, columns].
        
        :Parameters:
        
        group: int
            The group to be extracted.
        normalize: bool, optional
            If True, normalize the data so the average is 1.0.
            The default is False.
            
        :Returned:
        
        extracted_model: MiriDarkReferenceModel
            A MiriDarkReferenceModel with a single group extracted.
            
        """
        if self.data.ndim == 4:
        
            # Determine a new 3-D dark data model by extracting the
            # specified group.
            cube_data = self.data[:,group,:,:]
            cube_err = self.err[:,group,:,:]
            if normalize:
                posdata = np.where(cube_data > 0.0)
                factor = np.mean(cube_data[posdata])
#                 print("Normalizing by a factor of %f" % factor)
                if factor > 0.0:
                    cube_data = cube_data / factor
                    cube_err = cube_err / factor
            cube_dq = np.squeeze(self.dq)
            newmodel = MiriDarkReferenceModel( data=cube_data, err=cube_err,
                                               dq=cube_dq )
            # Copy the metadata to the new data model.
            # NOTE: The WCS metadata may be copied incorrectly.
            newmodel.copy_metadata( self )
            historystrg = "New DARK made by extracting group %d" % group
            if normalize:
                historystrg += " (normalized to 1.0)"
            newmodel.add_comment( historystrg )
            return newmodel
        
        elif self.data.ndim == 3:
        
            # Determine a new 2-D dark data model by extracting the
            # specified group.
            image_data = self.data[group,:,:]
            image_err = self.err[group,:,:]
            if normalize:
                posdata = np.where(image_data > 0.0)
                factor = np.mean(image_data[posdata])
#                 print("Normalizing by a factor of %f" % factor)
                if factor > 0.0:
                    image_data = image_data / factor
                    image_err = image_err / factor
            image_dq = np.squeeze(self.dq)
            newmodel = MiriDarkReferenceModel( data=image_data, err=image_err,
                                               dq=image_dq, averaged=True )
            # Copy the metadata to the new data model.
            # NOTE: The WCS metadata may be copied incorrectly.
            newmodel.copy_metadata( self )
            historystrg = "New DARK made by extracting group %d" % group
            if normalize:
                historystrg += " (normalized to 1.0)"
            newmodel.add_comment( historystrg )
            return newmodel
        
        else:
            raise TypeError("Data model has wrong number of dimensions")

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriDarkReferenceModel module.")

    PLOTTING = False
    SAVE_FILES = False

    data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    err3x3x2 = [err3x3, err3x3]
    err3x3x2x2 = [err3x3x2, err3x3x2]
    dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])
    dq3x3x2 = [dq3x3, dq3x3]
    dq3x3x2x2 = [dq3x3x2, dq3x3x2]
    data3x3x2 = [data3x3,data3x3]
    data3x3x2x2 = [data3x3x2,data3x3x2]

    print("Dark data with data + err + dq:")
    with MiriDarkReferenceModel(data=data3x3x2x2, err=err3x3x2x2, dq=dq3x3x2x2, \
                                dq_def=dark_reference_flags) \
            as testdata1:
        print(testdata1)
        if PLOTTING:
            testdata1.plot(description="testdata1")
        if SAVE_FILES:
            testdata1.save("test_darkreference_model1.fits", overwrite=True)
        del testdata1

    print("Dark data with data + err + dq:")
    with MiriDarkReferenceModel(data=data3x3x2x2, err=err3x3x2x2, dq=dq3x3x2x2,
                                dq_def=dark_reference_flags,
                                integration=42, fitted_after=39) as testdata2:
        print(testdata2)
        if PLOTTING:
            testdata2.plot(description="testdata2")
        if SAVE_FILES:
            testdata2.save("test_darkreference_model2.fits", overwrite=True)
        del testdata2

    print("Test finished.")
