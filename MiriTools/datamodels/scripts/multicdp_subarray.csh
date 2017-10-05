#!/bin/csh -f
#
# multicdp_subarray - Add missing subarray keywords to multiple CDPs.
#
# Specify a list of files, or wild cards, as script parameters.
#
# :Example:
#
# multicdp_subarray.csh *FLAT*.fits
#
# :History:
#
# 08 Dec 2015: Created.
#
# @author: Steven Beard (UKATC)
#
#===============================================================================
#
#
#  Apply the "cdp_add_subarray.py" command to each file in turn.
#
   foreach file ( $* )
      echo Adding subarray keywords to ${file} to make ${file}_subarray.fits ...
      cdp_add_subarray.py  ${file} ${file}_subarray.fits
   end
   echo Conversion finished.
