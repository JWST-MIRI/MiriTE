#!/bin/csh -f
#
# multicdp_remove_junk - Correct the BAND keyword within multiple CDPs.
#
# Specify a list of files, or wild cards, as script parameters.
#
# :Example:
#
# multicdp_remove_junk.csh *MIRIFU*.fits
#
# :History:
#
# 10 Aug 2018: Created.
#
# @author: Steven Beard (UKATC)
#
#===============================================================================
#
#
#  Apply the "cdp_remove_junk.py" command to each file in turn.
#
   foreach file ( $* )
      echo Removing junk from ${file} to make ${file}_nojunk.fits ...
      cdp_remove_junk.py  ${file} ${file}_nojunk.fits
   end
   echo Conversion finished.
