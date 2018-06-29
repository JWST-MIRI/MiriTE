#! /usr/bin/env python
"""

Simple examples of functions and classes in "dqflags" module.

:History:

25 Sep 2013: created (as test_dqflags.py)
03 Oct 2013: Converted to match changes to dqflags module and moved to "scripts".
29 Sep 2014: reserved_flags replaced by master_flags.
08 Sep 2015: Made compatible with Python 3.

@author: R. Azzollini (DIAS), Steven Beard (UKATC)

"""



#from pdb import set_trace as stop

import numpy as np


#from miri.dataproducts.dqflags import flags as fl
from miri.datamodels.dqflags import master_flags, FlagsTable

def mainx():
    
    flags_table = FlagsTable( master_flags )
    
    datashape = (10,15)  # shape of SCI and DQ arrays
    ndreg = [8,9,0,14]   # region where the CDP is not defined. 
                         # [x0,x1,y0,y1]
    
    satpix = (7,8)
    
    sci = np.ones(datashape,dtype='float32') # SCI array 
    dq = np.zeros(datashape,dtype='int32') # DQ array
    
    fksci = sci * 0.
    fksci[ndreg[0]:ndreg[1]+1,ndreg[2]:ndreg[3]+1] = np.nan
    ixndreg = np.where(np.isnan(fksci))
    
    # Nulling SCI where the CDP is not defined.
    sci[ixndreg] = 0. 
   
    # Saturated for a single pixel
    sci[satpix] = np.nan 
    
    # Setting "non defined" pixels as such (and discarding them)
    # Non Science data and Don't use it!
    dq[ixndreg] = flags_table.raiseflags(dq[ixndreg], ['NON_SCIENCE', 'DO_NOT_USE']) 
    
    # Setting "saturated" pixels as such (and discarding them)
    # Bad Solution and Don't use it        
    dq[np.where(np.isnan(sci))] = \
        flags_table.raiseflags(dq[np.where(np.isnan(sci))], ['SATURATED', 'DO_NOT_USE']) 
    
    print( 'pixels set as "saturated" are: ',\
        np.where(flags_table.test_flags_any(dq,'SATURATED')) )
    print( "Mean of pix's NOT 'Non Defined' = %.3f" % \
        sci[np.where(~flags_table.test_flags_any(dq, 'NON_SCIENCE'))].mean() )
    print( "Mean of pix's NOT 'Non Defined' AND NOT 'Saturated' = %.3f" %\
        (sci[np.where(~flags_table.test_flags_any(dq, 'NON_SCIENCE') &\
            ~flags_table.test_flags_any(dq, 'SATURATED'))]).mean() )
    print( "i.e. mean of pix's set as VALID = %.3f" % \
        (sci[np.where(~flags_table.test_flags_any(dq, 'DO_NOT_USE'))]).mean() )
    
    
    print( "Now we change 'saturated' pixel..." )
    
    sci[np.where(flags_table.test_flags_any(dq, 'SATURATED'))] = 1.
    dq[np.where(flags_table.test_flags_any(dq, 'SATURATED'))] = \
        flags_table.lowerflags(dq[np.where(flags_table.test_flags_any(dq,'SATURATED'))],\
            'SATURATED')
    
    print( 'pixel at (%i,%i) has now value %.f' % (satpix+(sci[satpix],)) )
    print( "Mean of pix's NOT 'Non Defined' = %.3f" % \
        sci[np.where(~flags_table.test_flags_any(dq,'NON_SCIENCE'))].mean() )
    

if __name__ == '__main__': mainx()
