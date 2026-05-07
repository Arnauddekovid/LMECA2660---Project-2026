from file_reader import *

M_05_static_until_001s = "M_0.5_static_until_0.01s.txt"
M_05_moving_until_001s = "M_0.5_moving_until_0.01s.txt"

M_05_static_until_01s = "M_0.5_static_until_0.1s.txt"

M_085_static_until_001s = "M_0.85_static_until_0.01s.txt"
M_085_moving_until_001s = "M_0.85_moving_until_0.01s.txt"

M_085_moving_until_01s = "M_0.85_moving_until_0.1s.txt"

M_085_static_until_01s = "M_0.85_static_until_0.1s.txt"
M_085_moving_until_01s = "M_0.85_moving_until_0.1s.txt"

growing_M_inf_until_01s = "M_0.5_growing_until_0.1s.txt"
growing_M_inf_until_05s = "M_0.5_growing_until_0.5s.txt"

test_static = "M_0.5_static_until_0.0001s.txt"
test_moving = "M_0.5_moving_until_0.001s.txt"

#_____

filename = M_085_static_until_001s

all_plots = False
if all_plots:
    Mach_number_over_time_animation(filename)
    Cp_over_time_animation(filename)
    normalized_total_enthalpy_over_time_animation(filename)
    normalized_total_pressure_over_time_animation(filename)
    normalized_entropy_over_time_animation(filename)
    normalized_vorticity_over_time_animation(filename)
else:
    #s_over_time_animation(filename)
    #s0_over_time_animation(filename)
    Mach_number_over_time_animation(filename)
    #evolution_animation_grid(filename)

