#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

This file contains all the global parameters describing the MIRI instrument,
its capabilities and modes of operation. It exists to ensure the 

09 Oct 2013: Created within datamodels.util
28 Jun 2018: Moved to miri.parameters
09 Aug 2018: Single and cross-dichroic pass-band names defined separately.
05 Sep 2018: FILTER, CHANNEL and BAND keywords are now compulsory in
             some CDP metadata.
17 Oct 2018: The string 'ANY' is no longer recommended within CDP metadata.
             'N/A' should be used instead.
12 Nov 2018: Added compulsory metadata for imager and MRS CDPs.
             Updated list of detector settings.
15 Nov 2018: Added CDP_USEAFTER_DICT.
13 Dec 2019: Removed obsolete FASTINTAVG readout mode.
05 Jun 2020: Reordered the subarray table to match the figure in PASP VIII.
28 Sep 2021: Added FASTR1, SLOWR1, and FASTGRPAVGn readout modes.
             Added new USEAFTER dates for specific subarrays.

@author: Steven Beard (UKATC)

"""

# --------------------------------------------------------------------------
# Dictionary of available detector readout modes.
#
# The tuple contains default values for:
# samplesum = number of A/D samples per reading
# sampleskip = number of A/D samples skipped before reading a pixel.
# refpixsampleskip = number of A/D samples skipped before reading a reference pixel.
# nframes = number of frames averaged per group.
# groupgap = number of frames dropped between groups.
# ngroups  = default number of readout groups per integration.
# nints    = default number of integrations per exposure.
# avggrps  = number of groups averaged to reduce data rate.
# avgints  = number of integrations averaged to reduce data rate.
#
#                                      r
#                                      e
#                                      f
#                                      p
#                                      i
#                                      x
#                                  s   s
#                              s   a   a
#                              a   m   m       g
#                              m   p   p   n   r   n       a   a
#                              p   l   l   f   o   g       v   v
#                              l   e   e   r   u   r   n   g   g
#                              e   s   s   a   p   o   i   g   i
#                              s   k   k   m   g   u   n   r   n
#                              u   i   i   e   a   p   t   p   t
#                              m   p   p   s   p   s   s   s   s
READOUT_MODE = {}
READOUT_MODE['SLOW'] =       (8,   1,  3,  1,  0, 10,  1,  1,  1)
READOUT_MODE['SLOWR1'] =     (8,   1,  3,  1,  0, 10,  1,  1,  1) # With extra reset
READOUT_MODE['FAST'] =       (1,   0,  3,  1,  0,  1, 10,  1,  1)
READOUT_MODE['FASTR1'] =     (1,   0,  3,  1,  0,  1, 10,  1,  1) # With extra reset
#READOUT_MODE['FASTINTAVG'] = (1,   0,  3,  1,  0,  1,  4,  1,  4)
READOUT_MODE['FASTGRPAVG']   = (1,   0,  3,  1,  0,  4,  1,  4,  1)
READOUT_MODE['FASTGRPAVG8']  = (1,   0,  3,  1,  0,  8,  1,  8,  1)
READOUT_MODE['FASTGRPAVG16'] = (1,   0,  3,  1,  0, 16,  1, 16,  1)
READOUT_MODE['FASTGRPAVG32'] = (1,   0,  3,  1,  0, 32,  1, 32,  1)
READOUT_MODE['FASTGRPAVG64'] = (1,   0,  3,  1,  0, 64,  1, 64,  1)
# The following four readout modes were used for MIRI testing only, and they
# will upset the JWST pipeline software if used.
# READOUT_MODE['SLOWINTAVG'] = (8,   1,  3,  1,  0,  1,  4,  1,  4)
# READOUT_MODE['SLOWGRPAVG'] = (8,   1,  3,  1,  0,  4,  1,  4,  1)
# READOUT_MODE['SLOWGRPGAP'] = (8,   1,  3,  4,  8,  4,  1,  1,  1)
# READOUT_MODE['FASTGRPGAP'] = (1,   0,  3,  4,  8,  4,  1,  1,  1)

# --------------------------------------------------------------------------
# Dictionary of available subarray options.
#
# REFERENCE:
# The Mid-Infrared Instrument for the James Webb Space Telescope,
# VIII: The MIRI Focal Plane System; M. E. Ressler et al.,
# Publications of the Astronomical Society of Pacific, Volume 127,
# Issue 953, pp. 675 (2015)
#
# Each tuple contains
# (first column, first row, number of columns, number of rows).
#                               s    s
#                               t    t     c
#                               a    a     o
#                               r    r     l
#                               t    t     u     r
#                               c    r     m     o
#                               o    o     n     w
#                               l    w     s     s
SUBARRAY = {}
SUBARRAY['FULL'] =          (   1,   1, 1032, 1024 )
SUBARRAY['GENERIC'] =       (   1,   1, 1032, 1024 )
SUBARRAY['BRIGHTSKY'] =     ( 457,  51,  512,  512 )
SUBARRAY['SUB256'] =        ( 413,  51,  256,  256 )
SUBARRAY['SUB128'] =        (   1, 889,  136,  128 )
SUBARRAY['SUB64'] =         (   1, 779,   72,   64 )
SUBARRAY['SLITLESSPRISM'] = (   1, 529,   72,  416 )
SUBARRAY['MASK1065'] =      (   1,  19,  288,  224 )
SUBARRAY['MASK1140'] =      (   1, 245,  288,  224 )
SUBARRAY['MASK1550'] =      (   1, 467,  288,  224 )
SUBARRAY['MASKLYOT'] =      (   1, 717,  320,  304 )

# --------------------------------------------------------------------------
# Lists of selections, as defined in
# http://jwst-reffiles.stsci.edu/source/required_keywords.html
#
# Which MIRI model have the data come from (verification model, JPL model or
# flight model)? 'VM' and 'JPL' will only be found in very old data.
MIRI_MODELS = ['VM', 'JPL', 'FM']

# The 3 detectors installed in the MIRI instrument.
MIRI_DETECTORS = ['MIRIMAGE', 'MIRIFULONG', 'MIRIFUSHORT']
MIRI_DETECTORS_EXTRAS = ['IM', 'LW', 'SW'] # For backwards compatibility only.

# Detector electronic settings used. This rarely changes.
MIRI_SETTINGS = ['RAL1', 'JPL1', 'JPL2', 'JPL3',
                 'RUN1', 'RUN2', 'RUN3', 'RUN4', 'RUN5', 'RUN6']

# Available MIRI readout modes.
MIRI_READPATTS = ['SLOW', 'SLOWR1', 'FAST', 'FASTR1', 'FASTGRPAVG',
                  'FASTGRPAVG8', 'FASTGRPAVG16', 'FASTGRPAVG32',
                  'FASTGRPAVG64']

# Available MIRI subarray modes. Compare with the SUBARRAY dictionary
# defined above. Also = list(SUBARRAY.keys()) - ['FULL', 'GENERIC']
MIRI_SUBARRAYS = ['BRIGHTSKY', 'SUB256', 'SUB128', 'SUB64', 'SLITLESSPRISM',
                  'MASK1065', 'MASK1140', 'MASK1550', 'MASKLYOT']

# Available MIRI MRS channels. Calibration data can be relevant for a single
# MRS channel or for a whole MRS detector containing two channels.
MIRI_CHANNELS_SINGLE = ['1', '2', '3', '4']
MIRI_CHANNELS_DOUBLE = ['12', '34']
MIRI_CHANNELS = MIRI_CHANNELS_SINGLE + MIRI_CHANNELS_DOUBLE

# Allowed settings for the two MIRI MRS dichroic grating wheels.
MIRI_DGAA = ['SHORT', 'MEDIUM', 'LONG']
MIRI_DGAB = ['SHORT', 'MEDIUM', 'LONG']

# Single passband settings when the dichroic wheels have the same setting.
MIRI_BANDS_SINGLE = ['SHORT', 'MEDIUM', 'LONG']
# Cross-dichroic settings when the dichroic wheels have different settings.
MIRI_BANDS_CROSS = ['SHORT-MEDIUM', 'SHORT-LONG',
                    'MEDIUM-SHORT', 'MEDIUM-LONG',
                    'LONG-SHORT',   'LONG-MEDIUM']
# The BAND keyword can contain any single or cross-dichroic passband.
MIRI_BANDS = MIRI_BANDS_SINGLE + MIRI_BANDS_CROSS
MIRI_BANDS_EXTRAS = ['A', 'B', 'C'] # For backwards compatibility only.

# Available MIRI imager filters.
# Note that 'FLENS' and 'F2550WR' are allowed values for the filter metadata,
# but there are no CDP files for these filters.
MIRI_FILTERS = ['F560W','F770W','F1000W','F1130W', 'F1280W','F1500W','F1800W',
                'F2100W', 'F2550W', 'F2550WR','F1065C', 'F1140C', 'F1550C',
                'F2300C','P750L','FLENS', 'FND', 'OPAQUE']

# --------------------------------------------------------------------------
# Rules for testing compulsory CDP metadata
CDP_METADATA = [['TELESCOP', 'JWST'],
                ['INSTRUME', 'MIRI'],
                ['MODELNAM', MIRI_MODELS + ['N/A']],
                ['DETECTOR', MIRI_DETECTORS + ['N/A']],
                ['DETSETNG', MIRI_SETTINGS + ['N/A']],
                ['READPATT', MIRI_READPATTS + ['N/A']],
                ['SUBARRAY', MIRI_SUBARRAYS + ['FULL', 'GENERIC', 'N/A']],
                ['FASTAXIS', 1],
                ['SLOWAXIS', 2],
                ['PEDIGREE', ['FLIGHT', 'GROUND', 'DUMMY', 'SIMULATION']],
                ['USEAFTER', []],  # Empty list means any value accepted.
                ['DESCRIP', []],  # Empty list means any value accepted.
                ['AUTHOR', []],  # Empty list means any value accepted.
                ['DATE', []],  # Empty list means any value accepted.
                ['VERSION', []],  # Empty list means any value accepted.
                ]
# The following keywords are compulsory for IMAGER and LRS CDPs
CDP_METADATA_IMAGER = [['FILTER', MIRI_FILTERS + ['N/A']]]
# The following keywords are compulsory for MRS CDPs
CDP_METADATA_MRS = [['CHANNEL', MIRI_CHANNELS + ['N/A']],
                ['BAND', MIRI_BANDS + ['N/A']],
                ]
CDP_METADATA_SUBSET = [['FILTER', MIRI_FILTERS + ['N/A']],
                ['CHANNEL', MIRI_CHANNELS + ['N/A']],
                ['BAND', MIRI_BANDS + ['N/A']],
                ]

# Additional compulsory metadata for non-GENERIC subarrays
CDP_SUBARRAY = [['SUBSTRT1', []], # Empty list means any value accepted.
                ['SUBSTRT2', []], # Empty list means any value accepted.
                ['SUBSIZE1', []], # Empty list means any value accepted.
                ['SUBSIZE2', []], # Empty list means any value accepted.
                ]
# Keywords used in HISTORY records
CDP_HISTORY = ['DOCUMENT', 'SOFTWARE', 'DATA USED', 'DIFFERENCES']

# A list of common dates for the reference file USEAFTER keyword
CDP_USEAFTER_DICT = {'DEFAULT'   : '2000-01-01T00:00:00',
                     'FM'        : '2011-05-01T00:00:00',
                     'CV1'       : '2013-09-26T05:00:00',
                     'CV2'       : '2013-09-26T05:00:00',
                     'CV3'       : '2013-09-26T05:00:00',
                     'CV3_BURST' : '2015-08-01T00:00:00',
                     'CDP1'      : '2013-03-01T00:00:00',
                     'SUB128'    : '2013-09-26T05:00:00',
                     'SUB64'     : '2013-09-26T05:00:00',
                     'SLITLESSPRISM' : '2013-09-26T05:00:00',
                     'MASK1065'  : '2013-09-26T05:00:00',
                     'MASK1140'  : '2013-09-26T05:00:00',
                     'MASK1550'  : '2013-09-26T05:00:00',
                     'MASKLYOT'  : '2013-09-26T05:00:00',
                     'CDP2'      : '2013-11-01T00:00:00',
                     'CDP3'      : '2014-11-01T00:00:00',
                     'CDP4'      : '2015-06-15T00:00:00',
                     'BRIGHTSKY' : '2015-08-01T00:00:00',
                     'SUB256'    : '2015-08-01T00:00:00',
                     'CDP5'      : '2015-12-18T00:00:00',
                     'CDP6'      : '2016-06-30T00:00:00',
                     'CDP7'      : '2018-12-01T00:00:00',
                     'CDP8'      : '2021-08-09T00:00:00',
                     'CDP9'      : '2022-01-01T00:00:00',
                    }

#
# Dictionary of the relationship between known detector settings
#                         Name -> (FRMRSETS, ROWRSETS, RPCDELAY)
DETECTOR_SETTINGS_DICT = {'RAL1': (0, 3, 24),
                          'JPL1': (3, 4, 90)}


if __name__ == '__main__':
    print("Testing MIRI parameters:")
    print(72 * "-")
    print("MIRI_MODELS=", MIRI_MODELS)
    print("MIRI_DETECTORS=", MIRI_DETECTORS)
    print("MIRI_SETTINGS=", MIRI_SETTINGS)
    print("MIRI_READPATTS=", MIRI_READPATTS)
    print("MIRI_SUBARRAYS=", MIRI_SUBARRAYS)
    print("MIRI_CHANNELS=", MIRI_CHANNELS)
    print("MIRI_BANDS=", MIRI_BANDS)
    print("MIRI_FILTERS=", MIRI_FILTERS)
    
    print(72 * "-")
    print("CDP_METADATA=")
    for item in CDP_METADATA:
        print("   %s" % item)
    print("CDP_SUBARRAY=", CDP_SUBARRAY)
    print("CDP_HISTORY=", CDP_HISTORY)

    print(72 * "-")
    print("Test finished.")
