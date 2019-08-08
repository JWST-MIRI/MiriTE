#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
import apt_parser
import logging
apt_parser.init_log()

LOG = logging.getLogger('parse_apt')

folder = "../../xml"
files = glob.glob(os.path.join(folder, "*.aptx"))

(APT_sim, mem_volume, time_volume) = apt_parser.analyse_apt_list(files)

print("mem_volume")
print(mem_volume)
print("time_volume")
print(time_volume)

# Without parallel
all_time = 0
all_max_ram = 0
for car in time_volume.keys():
    total_time = time_volume[car].sum()
    max_ram = mem_volume[car].max()

    all_time += total_time
    all_max_ram = max(all_max_ram, max_ram)

    print("{} ({}.aptx): Total Time = {:.1f} h ; Max memory = {:.1f} GB".format(car, apt_parser.constants.MIRI_STSCI_ID[car], total_time, max_ram))

days = int(all_time // 24)
hours = all_time - 24 * days
print("Total time for all files: {:d} d {:.1f} h ; max RAM needed: {:.1f} GB".format(days, hours, all_max_ram))
# JwstProposal/DataRequests/ObservationGroup