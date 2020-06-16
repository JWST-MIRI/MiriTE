#!/bin/sh
#
python scasim.py --detector MIRIFUSHORT --rdmode $1 --ngroups 100 --nints 1 --crmode NONE --norefpixels --nobadpixels --nodark --noflat --nolinearity --nodrifts --nolatency --scale 20.0 ConstantIllumination.fits $2 $3 $4
