#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Running all simulations related to a specific APT file
"""

import logging
import apt_parser
import apt_parser.mirisim

LOG = logging.getLogger('main')

# JWST imports
from mirisim import skysim
from mirisim import config_parser as c  # SimConfig, SimulatorConfig, SceneConfig

# To hide MIRISim INFO level
extra_config = {"loggers":
            {

                "mirisim":
                {
                        "level": "WARNING",
                    },
            },}

apt_parser.init_log(log="mirisim.log", stdout_loglevel="INFO", file_loglevel="DEBUG", extra_config=extra_config)

# We prepare the scene we want to model
background = skysim.Background(level='low', gradient=0., pa=0.0, centreFOV=(0., 0.))

annulus = skysim.Skycube('SN1987A.fits', center='yes')

star_sed = skysim.BBSed(Temp=90., wref=10.6, flux=200.)
star = skysim.Point(Cen=(0.,0.))
star.set_SED(star_sed)

# Coordinates are extracted from visible_w3_r_2018-07-08_drz_al2.fits, then the relative coordinates
# are computed by the file get_coordinates.py
# Fluxes are:
# S2 : B2III, Vmag = 15, T = 22000 K
# S3 : A8V, Vmag = 16, T = 7500 K
# Sx : K5V, Vmag = 17.7, T = 4330 K
sed_S2 = skysim.BBSed(Temp=22000., wref=0.55, flux=3.78e3)
S2 = skysim.Point(Cen=(-1.896, 2.140))
S2.set_SED(sed_S2)

sed_S3 = skysim.BBSed(Temp=7500., wref=0.55, flux=1.51e3)
S3 = skysim.Point(Cen=(1.487, -0.824))
S3.set_SED(sed_S3)

sed_SX = skysim.BBSed(Temp=4330., wref=0.55, flux=0.31e3)
SX = skysim.Point(Cen=(-1.869, -1.004))
SX.set_SED(sed_SX)

targets = [annulus, star, S2, S3, SX]

# scene.ini
scene_config = c.SceneConfig.makeScene(loglevel=1, background=background, targets=targets)

observations = apt_parser.parse_apt("1232.aptx")

apt_parser.mirisim.run(observations, scene_config)

