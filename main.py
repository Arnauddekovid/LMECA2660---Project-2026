from file_reader import *

filename_M050_t01 = "M_inf = 0.500000 until t_end = 0.100000.txt"
filename_M050_t001 = "M_inf = 0.500000 until t_end = 0.010000.txt"
filename_M050_t0001 = "M_inf = 0.500000 until t_end = 0.001000.txt"
filename_M050_t00001 = "M_inf = 0.500000 until t_end = 0.000100.txt"
filename_M050_t000001 = "M_inf = 0.500000 until t_end = 0.000010.txt"


filename_M085_t1 = "M_inf = 0.850000 until t_end = 1.000000.txt"
filename_M085_t01 = "M_inf = 0.850000 until t_end = 0.100000.txt"
filename_M085_t001 = "M_inf = 0.850000 until t_end = 0.010000.txt"
filename_M085_t0001 = "M_inf = 0.850000 until t_end = 0.001000.txt"
filename_M085_t00001 = "M_inf = 0.850000 until t_end = 0.000100.txt"
filename_M085_t000001 = "M_inf = 0.850000 until t_end = 0.000010.txt"


data_D = "D  M_inf = 0.500000 until t_end = 0.000100.txt"

delete_txt_file_after_animation = False
#plot_grid(filename_M050_t000001)

evolution_animation(filename_M085_t01, delete_txt_file_after_animation)

print("finished")

#plot_other_data(data_D, "D ")
