# -*- coding: utf-8 -*-

"""

This module contains utilities for testing the CDP fetching utilites.

:Reference:

http://miri.ster.kuleuven.be/bin/view/Internal/CalDataProducts

:History:

31 Oct 2016: Created

Steven Beard (UKATC), Vincent Geers (UKATC)

"""

# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Python logging facility.
import logging
logging.basicConfig(level=logging.INFO)   # Choose ERROR, WARN, INFO or DEBUG 
LOGGER = logging.getLogger("miri.cdp_tests") # Get a default parent logger

from miri.datamodels.cdplib import get_cdp as miri_get_cdp

#
# A minimal test and some examples of how to use the above utilities
# are run when this file is executed as a main program.
#
if __name__ == '__main__':
    VERBOSE = True
    PLOTTING = False
    
    import time    
    print("\nTesting use of get_cdp")

    dm1 = miri_get_cdp('FLAT', detector='MIRIFULONG', band='SHORT',
                       ftp_path='CDPSIM', ftp_user='miriuser',
                       ftp_passwd='jWsTm1R1', logger=LOGGER)
    print( dm1 )
    dm2 = miri_get_cdp('FLAT', detector='MIRIFULONG', band='MEDIUM',
                       ftp_path='CDPSIM', ftp_user='miriuser',
                       ftp_passwd='jWsTm1R1', logger=LOGGER)
    print ( dm2 )
    dm3 = miri_get_cdp('FLAT', detector='MIRIFULONG', band='LONG',
                       ftp_path='CDPSIM', ftp_user='miriuser',
                       ftp_passwd='jWsTm1R1', logger=LOGGER)
    print( dm3 )

    #The first two will fetch 34SHORT and 34MEDIUM, as I expected, but the
    # third fetches 34SHORT instead of 34LONG for me.

    #Explicitly setting channel will circumvent it:
    dm4 = miri_get_cdp('FLAT', detector='MIRIFULONG', channel='34', band='LONG',
                       ftp_path='CDPSIM', ftp_user='miriuser',
                       ftp_passwd='jWsTm1R1')
    print( dm4 )
    
    #Is it supported to call get_cdp without "channel" for MRS?

    print("Test finished.")
