#!/bin/csh -f
#
# multicdp_band - Correct the BAND keyword within multiple CDPs.
#
# Specify a list of files, or wild cards, as script parameters.
#
# :Example:
#
# multicdp_band.csh *MIRIFU*.fits
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
#  Apply the "cdp_correct_band.py" command to each file in turn.
#
   foreach file ( $* )
      echo Correcting band keyword in ${file} to make ${file}_newband.fits ...
      cdp_correct_band.py  ${file} ${file}_newband.fits
   end
   echo Conversion finished.
