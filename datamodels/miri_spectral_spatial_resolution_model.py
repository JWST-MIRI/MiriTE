#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines the MIRI
spectral and spatial resolution models.

:Reference:

The STScI jwst_lib documentation.

:History:

30 Aug 2016: Created.
12 Oct 2016: Renamed class from MiriMrsSpectralResolutionModel to
             MiriMrsResolutionModel and data type from SPECRES to
             RESOL to reflect the fact that this data model
             describes both the spectral and spatial resolution of
             the MRS. psf_alpha extension renamed to psf_fwhm_alpha
             and psf_beta extension renamed to psf_fwhm_beta.
17 Oct 2016: Data model schema modified to match above changes.
             SUB_BAND column added to RESOLVING_POWER table.
             Identical copies of the same CDP file generated for
             MIRIFULONG and MIRIFUSHORT.
06 Dec 2016: SUB_BAND column added to YAML version of schema.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.

@author: Steven Beard (UKATC), Alvaro Labiano (CSIC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

#import warnings
import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriMrsResolutionModel']


class MiriMrsResolutionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI spectral resolution table,
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
        
    resolving_power: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the spectral
        resolving power table;
        Or: A numpy record array containing the same information as above.
        A spectral resolution table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    psf_fwhm_alpha: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the spatial
        resolution alpha table;
        Or: A numpy record array containing the same information as above.
        A spatial resolution table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    psf_fwhm_beta: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the spatial
        resolution beta table;
        Or: A numpy record array containing the same information as above.
        A spatial resolution table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst_lib documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_spectral_spatial_resolution_mrs.schema.yaml"
    fieldnames_resolving = ('SUB_BAND', 'R_CENTRE', 'R_A_LOW', 'R_B_LOW', 'R_C_LOW',
                            'R_A_HIGH', 'R_B_HIGH', 'R_C_HIGH',
                            'R_A_AVG', 'R_B_AVG', 'R_C_AVG')
    fieldnames_alpha = ('A_CUTOFF','A_A_SHORT', 'A_B_SHORT',
                        'A_A_LONG', 'A_B_LONG')
    fieldnames_beta = ('B_CUTOFF', 'B_A_SHORT', 'B_B_SHORT',
                       'B_A_LONG', 'B_B_LONG')
    
    def __init__(self, init=None, resolving_power=None, psf_fwhm_alpha=None,
                 psf_fwhm_beta=None, **kwargs):
        """
        
        Initialises the MiriMrsResolutionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsResolutionModel, self).__init__(init=init, **kwargs)

        # Data type is spectral resolution.
        self.meta.model_type = 'RESOL'
        self.meta.reftype = 'RESOL'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'
        
        if resolving_power is not None:
            try:
                self.resolving_power = resolving_power
            except (ValueError, TypeError) as e:
                strg = "resolving_power must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if psf_fwhm_alpha is not None:
            try:
                self.psf_fwhm_alpha = psf_fwhm_alpha
            except (ValueError, TypeError) as e:
                strg = "psf_fwhm_alpha must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if psf_fwhm_beta is not None:
            try:
                self.psf_fwhm_beta = psf_fwhm_beta
            except (ValueError, TypeError) as e:
                strg = "psf_fwhm_beta must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        
#         # Copy the table column units, if defined.
#         resolving_power_units = self.set_table_units('resolving_power')
#         psf_fwhm_alpha_units = self.set_table_units('psf_fwhm_alpha')
#         psf_fwhm_beta_units = self.set_table_units('psf_fwhm_beta')
        
    def __str__(self):
        """
        
        Display the contents of the spectral resolution object
        as a readable string.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')

        # Describe the spectral resolution tables
        if self.resolving_power is not None:
            strg += self.get_data_str('resolving_power', underline=True, underchar="-")
        if self.psf_fwhm_alpha is not None:
            strg += self.get_data_str('psf_fwhm_alpha', underline=True, underchar="-")
        if self.psf_fwhm_beta is not None:
            strg += self.get_data_str('psf_fwhm_beta', underline=True, underchar="-")
        return strg


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriMrsResolutionModel module.")
    
    PLOTTING = False
    SAVE_FILES = True

    # Sub-bands
    sb = ['1SHORT', '1MEDIUM', '1LONG',
          '2SHORT', '2MEDIUM', '2LONG',
          '3SHORT', '3MEDIUM', '3LONG',
          '4SHORT', '4MEDIUM', '4LONG'
          ]

    # Centros
    cntr = [5.3, 6.15, 7.1, 8.15, 9.4, 10.9, 12.5, 14.5, 16.75, 19.2, 22.25, 26.05]
    
    # Low limits
    al = [2745.211671, 2737.58174, 2691.826643, 2716.566802, 2575.145064, 2563.664138, 2469.622611, 1864.309643, 2071.315475, 1899.989987, 1547.295843, 1220.329908]
    bl = [541.3274075, 506.8427022, 393.1504096, 324.8195469, 280.4117705, 257.0795746, 128.1952696, -24.10526842, 268.5610664, -50.94618217, 12.23096891, -146.1678896]
    cl = [310.8853084, 250.7537815, 35.05343649, -259.5669385, -182.4385113, -122.1030201, -279.8033748, -9.54407769, -23.40212084, -81.27269184, -12.23916046, -60.13030606]
    
    # High limits
    ah = [4063.403427, 3934.231411, 3742.01575, 4748.404518, 4227.307398, 4067.791865, 3642.877564, 2768.145139, 2641.910569, 2158.932162, np.nan, np.nan]
    bh = [193.251713, 151.0248429, 79.69070813, 476.6602166, 420.9441793, 227.7079654, 326.3927585, -86.76805885, 363.0089443, -70.19831027, np.nan, np.nan]
    ch = [-191.9995708, -418.0899389, -360.9548965, -1804.289718, -611.3820685, -569.5866384, -601.634922, -206.1799239, -45.5071111, -109.3828766, np.nan, np.nan]
    
    # Use values
    au = [3404.307549, 3335.906575, 3216.921197, 3732.48566, 3401.226231, 3315.728002, 3056.250087, 2316.227391, 2356.613022, 2029.461074, 1547.295843, 1220.329908]
    bu = [367.2895603, 328.9337726, 236.4205589, 400.7398818, 350.6779749, 242.39377, 227.2940141, -55.43666364, 315.7850053, -60.57224622, 12.23096891, -146.1678896]
    cu = [59.44286876, -83.66807868, -162.95073, -1031.928328, -396.9102899, -345.8448292, -440.7191484, -107.8620008, -34.45461597, -95.32778422, -12.23916046, -60.13030606]

    resolving_power = []
    for ii in range(0, len(cntr)):
        row = (sb [ii], cntr[ii], al[ii], bl[ii], cl[ii], ah[ii], bh[ii], ch[ii],
               au[ii], bu[ii], cu[ii])
        resolving_power.append(row)

##    An alternative way of defining the same information.
#     resolving_power = [( '1SHORT', 5.3, 2745.211671, 541.3274075, 310.8853084,
#                               4063.403427, 193.251713, -191.9995708, 
#                               3404.307549, 367.2895603, 59.44286876),
#                        ('1MEDIUM', 6.15, 2737.58174, 506.8427022, 250.7537815,
#                               3934.231411, 151.0248429, -418.0899389,
#                               3335.906575, 328.9337726, -83.66807868)]

    psf_fwhm_alpha = [(8.0, 0.31, 0.0, 0.0, 0.03875)]
    psf_fwhm_beta = [(8.0, 0.0, 0.03875, 0.0, 0.03875)]

    print("\nMRS spectral and spatial resolution model for MIRIFUSHORT:")
    with MiriMrsResolutionModel( resolving_power=resolving_power,
                                        psf_fwhm_alpha=psf_fwhm_alpha,
                                        psf_fwhm_beta=psf_fwhm_beta ) \
                                    as testspecres1:
        testspecres1.set_referencefile_metadata( author='Alvaro Labiano',
                    pedigree='GROUND', useafter='DEFAULT',
                    description='MIRI MRS Spectral and Spatial Resolution CDP')
        testspecres1.add_referencefile_history(
                    document='MIRI-TN-00005-ETH Issue 1.0',
                    software='IDL and Python',
                    dataused='Derived from FM data',
                    differences='N/A')
        testspecres1.set_instrument_metadata('MIRIFUSHORT', modelnam='FM',
                    detsetng='ANY', filt='ANY', channel='12', band='ANY')
        testspecres1.set_subarray_metadata('GENERIC')
        testspecres1.meta.exposure.readpatt = 'ANY'
        testspecres1.meta.exposure.type = 'MIR_MRS'
        print(testspecres1)
        if PLOTTING:
            testspecres1.plot(description="testspecres1")
        if SAVE_FILES:
            testspecres1.save("MIRI_FM_MIRIFUSHORT_12_RESOL_06.00.00.fits",
                              overwrite=True)
        del testspecres1

    print("\nMRS spectral and spatial resolution model for MIRIFULONG:")
    with MiriMrsResolutionModel( resolving_power=resolving_power,
                                        psf_fwhm_alpha=psf_fwhm_alpha,
                                        psf_fwhm_beta=psf_fwhm_beta ) \
                                    as testspecres2:
        testspecres2.set_referencefile_metadata( author='Alvaro Labiano',
                    pedigree='GROUND', useafter='DEFAULT',
                    description='MIRI MRS Spectral and Spatial Resolution CDP')
        testspecres2.add_referencefile_history(
                    document='MIRI-TN-00005-ETH Issue 1.0',
                    software='IDL and Python',
                    dataused='Derived from FM data',
                    differences='N/A')
        testspecres2.set_instrument_metadata('MIRIFULONG', modelnam='FM',
                    detsetng='ANY', filt='ANY', channel='34', band='ANY')
        testspecres2.set_subarray_metadata('GENERIC')
        testspecres2.meta.exposure.readpatt = 'ANY'
        testspecres2.meta.exposure.type = 'MIR_MRS'
        print(testspecres2)
        if PLOTTING:
            testspecres2.plot(description="testspecres2")
        if SAVE_FILES:
            testspecres2.save("MIRI_FM_MIRIFULONG_34_RESOL_06.00.00.fits",
                              overwrite=True)
        del testspecres2
        
    print("Test finished.")
