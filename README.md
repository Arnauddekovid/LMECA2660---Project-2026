========================================

LMECA2660 - PROJET 2026 : Arnaud Saussez

========================================

Prerequisites :

1 - Having some storage or use a SSD ( and time): 
    we store the data inside the binaryFiles folder, for some huge computations with high resolution, we can get a 25Go file, you can use an alternative folder by changing the storage_hdd value in the main.c file. The location of this alt folder can be changed inside the generate_filename() function at line 58 inside the file_writer.h file.

========================================

Usage :

2 - Go to the main.c file and use the classic cases provided inside the different_cases.h file:

    static_airfoil : the airfoil is static and can have an angle of attack
    oscillating_airfoil : the airfoil is moving at a frequency and rise to a certain value of amplitude
    rising_mach : static airfoil but the mach number is rising, useful to simulate acceleration and see the formation of shocks
                    you must start with a Mach_init close to 1 to avoid the non implemeted bounds switch 

    run the code with the following line to be way faster than the non optimised executor
            gcc main.c -O3 -march=native -o exec     

3 - plotting results

    In the main.py file, you can plot several plots but its better to use the all_plots_animations(), which is very faster than using every plotting function independenlty. the all plots animation provides you the 33fps videos of:

        - the Mach time evolution
        - the pressure Coefficient Cp time evolution
        - the normalized total enthalpy time evolution
        - the normalized total pressure time evolution
        - the normalized enthalpy time evolution
        - the normalized vorticity time evolution
        - the Mach and Cp on walls time evolution
        - the lift and drag coefficients time evolution
    
    
