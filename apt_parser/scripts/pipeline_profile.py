
import numpy as np

time = np.array([1408, 272.5, 132.8, 9.227, 5.579, 4.941, 2.502, 2.498, 0.5486, 0.3323, 0.278, 0.1141])

step = ["ramp_fit_step", "jump_step", "rscd_step", "ipc_step", "refpix_step", "linearity_step", "dark_current_step", "saturation_step", "dq_init_step", "lastframe_step", "firstframe_step", "gain_scale_step"]

total = time.sum()
total_time = 1850

percentage = time * 100 / total_time

# for (p, s) in zip(percentage, step):
#     print("{}: {:.2f}%".format(s, p))

print("\n".join(step))

for p in percentage:
    print("{:.2f}%".format(p))
