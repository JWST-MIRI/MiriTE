#!/bin/csh -f
#
# multicdp_filter_band - Add missing FILTER/CHANNEL/BAND within multiple CDPs.
#
# Specify a list of files, or wild cards, as script parameters.
#
# :Example:
#
# multicdp_filter_band.csh *MIRI*FLAT*.fits
#
# :History:
#
# 06 Sep 2018: Created.
#
# @author: Steven Beard (UKATC)
#
#===============================================================================
#
#
#  Apply the "cdp_add_filter_band.py" command to each file in turn.
#
   foreach file ( $* )
      echo Correcting filter/channel/band keywords in ${file} to make ${file}_newfilterband.fits ...
      cdp_add_filter_band.py  ${file} ${file}_newfilterband.fits
   end
   echo Conversion finished.
