//
//  creates different_cases
//  MECA2660 projet 2026


#ifndef different_cases_h
#define different_cases_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  
#include <complex.h>

#include "grid_and_data.h"
#include "file_writer.h"
#include "solver.h"
#include "debug_functions.h"



// static airfoil

void static_arfoil(double M_inf,double t_end, double angle_of_attack, double number_of_written_steps, double resolution, bool hdd_storage){
    /*
        pre :
            M_inf : the value of the Mach Numer at the infinity
            t_end : the time we stop the computations
            angle_of_attack : angle of attack of the airfoil, must be in degrees
            number_of_written_steps : how many steps we save on the file, 100 is equavelent of 3 seconds of video, 1000 of 33 
            resolution : 0.0 for low resolution, 0.5 for mid and 1.0 for high res. the high resolution is super heavu (25 Go for 1000number_of_written_steps )
            hdd_storage : true if the files will be too heavy and stores it on the hdd, the location of the hdd mus tbe change in the generate_filename() function in the file writer file
        post :
            geenrates a file of the problem, with a static airfoil, with some angle of attack. 
    */
    double grow_rate_M_inf = 0.0;

    double amplitude_oscillation = 0.0;
    double frequency_oscillation = 0.;


    solver_all_steps(M_inf, grow_rate_M_inf, t_end, number_of_written_steps,resolution,  amplitude_oscillation, frequency_oscillation,angle_of_attack, hdd_storage);
}

void moving_airfoil(double M_inf, double t_end, double amplitude_oscillation, double frequency_oscillation,  double number_of_written_steps, double resolution, bool hdd_storage){
    /*
        pre :
            M_inf : the value of the Mach Numer at the infinity
            t_end : the time we stop the computations
            amplitude_oscillation : amplitude of the oscillation, the higheest value of the leading/trailing edges
            frequency_oscillation : the frequency of the oscillations 
            number_of_written_steps : how many steps we save on the file, 100 is equavelent of 3 seconds of video, 1000 of 33 
            resolution : 0.0 for low resolution, 0.5 for mid and 1.0 for high res. the high resolution is super heavu (25 Go for 1000number_of_written_steps )
            hdd_storage : true if the files will be too heavy and stores it on the hdd, the location of the hdd mus tbe change in the generate_filename() function in the file writer file
        post :
            geenrates a file of the problem, with a moving airfoil, with some frequency of oscillation
    */
    double grow_rate_M_inf = 0.0;
    double angle_of_attack = 0.; // en degré

    solver_all_steps(M_inf, grow_rate_M_inf, t_end, number_of_written_steps,resolution,  amplitude_oscillation, frequency_oscillation,angle_of_attack, hdd_storage);

}

void rising_mach(double M_inf_init, double M_inf_end, double t_end,  double angle_of_attack, double number_of_written_steps, double resolution, bool hdd_storage){
    /* NOT WORKING NOW
        pre :
            M_inf : the value of the Mach Numer at the infinity at the beginning of the simulation
            M_inf_end : the value of the Mach Numer at the infinity at the end of the simulation
            t_end : the time we stop the computations
            angle_of_attack : angle of attack of the airfoil, must be in degrees
            number_of_written_steps : how many steps we save on the file, 100 is equavelent of 3 seconds of video, 1000 of 33 
            resolution : 0.0 for low resolution, 0.5 for mid and 1.0 for high res. the high resolution is super heavu (25 Go for 1000number_of_written_steps )
            hdd_storage : true if the files will be too heavy and stores it on the hdd, the location of the hdd mus tbe change in the generate_filename() function in the file writer file
        post :
            geenrates a file of the problem, with a  with a static airfoil, with some angle of attack. 
    */
    double grow_rate_M_inf = (M_inf_end-M_inf_init)/t_end;

    double amplitude_oscillation = 0.0;
    double frequency_oscillation = 0.;
    //bool hdd_storage = false; // true pour stocker les données sur le disque dur et pas sur un ordi, changer la destination du disque dure dans file_writer

    solver_all_steps(M_inf_init, grow_rate_M_inf, t_end, number_of_written_steps,resolution,  amplitude_oscillation, frequency_oscillation,angle_of_attack, hdd_storage);

}


#endif 