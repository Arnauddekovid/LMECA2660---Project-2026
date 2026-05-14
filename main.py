from file_reader import *


filename = "binaryFiles/high_res_M_inf_0.85_oscillating_mu_0.05_until_0.1s.bin"

all_plots = True
if all_plots:
    all_plots_animations(filename)
else:
    normalized_vorticity_over_time_animation(filename)
