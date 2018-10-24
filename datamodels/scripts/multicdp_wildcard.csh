#!/bin/csh -f
#
# multicdp_wildcard - Correct the BAND keyword within multiple CDPs.
#
# Specify a list of files, or wild cards, as script parameters.
#
# :Example:
#
# multicdp_wildcard.csh *MIRIFU*.fits
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
#  Apply the "cdp_correct_wildcard.py" command to each file in turn.
#
   foreach file ( $* )
      echo Correcting wildcard keyword in ${file} to make ${file}_wildcard.fits ...
      cdp_correct_wildcard.py  ${file} ${file}_wildcard.fits
   end
   echo Conversion finished.
