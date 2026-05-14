#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  
#include <complex.h>
#include <time.h>

#include "different_cases.h"
#include "solver.h"


// gcc main.c -o exec                                                           550 secondes pour 1000 steps
// gcc main.c -O3 -march=native -o exec                                     120 secondes pour 1000 steps -> best

// gcc main.c -o exec -g -fsanitize=address -lm

int main(){

    bool hdd_storage = false; // if the data could be too heavy
    double resolution = 1.; // 0.0 for a low resolution, 0.5 for a midgle one, 1. for high resolution
    double t_end =  .25;
    double number_of_written_steps = 1000. ; // 33 secondes de vidéo pour 1000 steps, avec la grosse grid,  fait environ 25Go, 3 secondes pour 100 stepsfait environ 2.5Go avec 0.01 en te_end

    double angle_of_attack = 0.; // must be in degrees
    double M_inf = .85;

    //static_arfoil(M_inf, t_end, angle_of_attack, number_of_written_steps, resolution, hdd_storage);

    double amplitude_oscillation = 0.1; // maximum value of the airfoil while oscillating
    t_end = .1;
    double frequency_oscillation = 10.*M_PI/t_end;
    moving_airfoil(M_inf, t_end, amplitude_oscillation, frequency_oscillation, number_of_written_steps, resolution, hdd_storage);

    t_end =  .1;
    number_of_written_steps = 1000. ; // 33 secondes de vidéo pour 1000 steps, avec la grosse grid,  fait environ 25Go, 3 secondes pour 100 stepsfait environ 2.5Go avec 0.01 en te_end
    resolution = .5; // 0.0 for a low resolution, 0.5 for a midgle one, 1. for high resolution
    double M_inf_init = .85;
    double M_inf_end = 3.;
    //rising_mach(M_inf_init , M_inf_end, t_end, angle_of_attack, number_of_written_steps, resolution, hdd_storage);



    return 0;
}
