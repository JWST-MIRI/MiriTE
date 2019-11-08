#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import glob
import os
from miri import apt_parser
from pipeline_parallel.fake_server import FakeServer
import logging
apt_parser.init_log()
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

LOG = logging.getLogger('parse_apt')


folder = "apt"
files = glob.glob(os.path.join(folder, "*.aptx"))
(APT_sim, mem_volume, time_volume) = apt_parser.analyse_apt_list(files)

# for car in mem_volume.keys():
#     runtime = predict_parallel_time(mem_volume[car], time_volume[car], total_ram=TOTAL_RAM, total_cpu=TOTAL_CPU)
#     print("With {} CPU and {} RAM, {} takes {:.1f} hours".format(TOTAL_CPU, TOTAL_RAM, car, runtime))


def optimize_server(rams_input, times_input, ylabel, title=None):
    server_cpus = np.arange(1,50)
    server_rams = [16, 32, 64, 128, 256, 512, 1024]  # GB
    runtimes = {}

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    for ram in server_rams:
        runtime = []

        for cpu in server_cpus:
            s = FakeServer(rams_input, times_input, total_ram=ram, total_cpu=cpu)
            r = s.predict_parallel_time()
            runtime.append(r)

        # NaN if not enough ram
        if np.isfinite(r).all():
            runtimes[ram] = runtime


    fig, ax = plt.subplots()
    runtime_keys = list(runtimes.keys())
    runtime_keys.sort()
    for i, ram in enumerate(server_rams):
        if ram in runtimes.keys():
            ax.semilogy(server_cpus, runtimes[ram], color=colors[i], label="{} GB".format(ram))
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Number of CPU")
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, which='major', color='#aaaaaa', linestyle='--')
    ax.yaxis.grid(True, which='major', color='#aaaaaa', linestyle='--')
    ax.yaxis.grid(True, which='minor', color='#cccccc', linestyle=':')
    ax.legend()
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=1)
    myyfmt = FormatStrFormatter("%.3g")
    ax.yaxis.set_major_formatter(myyfmt)
    ax.yaxis.set_minor_formatter(myyfmt)
    ax.tick_params(axis='both', which='minor', labelsize=7)

    if title:
        fig.suptitle(title)

    return fig


car = "MIRI-005"
miri_05_runtime = optimize_server(mem_volume[car], time_volume[car], ylabel="Total runtime [h]")
miri_05_runtime.savefig("server_time_{}.pdf".format(car), bbox_inches="tight", pad_inches=0.1)

# car = "MIRI-050"
# miri_50_runtime = optimize_server(mem_volume[car], time_volume[car], title="{}: Lots of medium files".format(car), ylabel="Total runtime [h]")
# miri_50_runtime.savefig("server_time_{}.pdf".format(car), bbox_inches="tight", pad_inches=0.1)

car = "MIRI-011"
miri_11_runtime = optimize_server(mem_volume[car], time_volume[car], ylabel="Total runtime [h]")
miri_11_runtime.savefig("server_time_{}.pdf".format(car), bbox_inches="tight", pad_inches=0.1)


plt.show()