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



void test_static(){    
    double M_inf = 0.5;
    double t_end =  1e-4;
    int step_max = (int) 1. ;

    double grow_rate_M_inf =0.; // on veut que M_inf augmente de  1 en t_end secondes, donc on a besoin d'un taux de croissance de 2/t_end jsp pq
    
    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.0; // poin
    double frequence_oscillation =0;

    struct grid_and_data grid = init_grid_and_data(M_inf, grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("static test done");
}

void test_moving(){    
    double M_inf = 0.5;
    double t_end =  1e-3;
    int step_max = (int) 1. ;

    double grow_rate_M_inf =0.; // on veut que M_inf augmente de  1 en t_end secondes, donc on a besoin d'un taux de croissance de 2/t_end jsp pq
    
    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.2; // poin
    double frequence_oscillation = M_PI_2/t_end;

    struct grid_and_data grid = init_grid_and_data(M_inf, grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("moving test done");

}


void growing_M_inf_until_01s(){
    double M_inf = 0.5;
    double t_end =  .1;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf =6./t_end; // on veut que M_inf augmente de  1 en t_end secondes, donc on a besoin d'un taux de croissance de 2/t_end jsp pq
    
    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.;
    double frequence_oscillation = 0.;
    
    struct grid_and_data grid = init_grid_and_data(M_inf, grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    grow_L_to_7(&grid);
    solver_all_steps(&grid);

    free_data(&grid);
    printf("growing_M_inf_until_01s done");

}

void growing_M_inf_until_05s(){
    double M_inf = 0.5;
    double t_end =  .5;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 4./t_end; // on veut que M_inf augmente de  1 en t_end secondes, donc on a besoin d'un taux de croissance de 2/t_end jsp pq
    
    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.;
    double frequence_oscillation = 0.;
    
    struct grid_and_data grid = init_grid_and_data(M_inf, grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    grow_L_to_7(&grid);
    solver_all_steps(&grid);

    free_data(&grid);
        
    printf("growing_M_inf_until_05s done");

}


// M_inf = 0.5


void M_05_static_until_001s(){
    double M_inf = 0.5;
    double t_end =  0.01;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave t_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.0;
    double frequence_oscillation = M_PI_2/t_end;

    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation, angle_of_attack);

    solver_all_steps(&grid);
    free_data(&grid);
    printf("M_05_moving_until_001s done");
}

void M_05_moving_until_001s(){
    double M_inf = 0.5;
    double t_end =  0.01;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave t_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.05;
    double frequence_oscillation = M_PI_2/t_end;

    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation, angle_of_attack);

    solver_all_steps(&grid);
    free_data(&grid);
    printf("M_05_moving_until_001s done");
}


void M_05_static_until_01s(){
    double M_inf = 0.5;
    double t_end =  0.1;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave t_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.;
    double frequence_oscillation = 0.;

    struct grid_and_data grid = init_grid_and_data(M_inf, grow_rate_M_inf,t_end, step_max, amplitude_oscillation, frequence_oscillation, angle_of_attack);

    solver_all_steps(&grid);
    free_data(&grid);
    printf("M_05_static_until_01s done");
}

void M_05_moving_until_01s(){
    double M_inf = 0.5;
    double t_end =  0.1;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave t_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.05;
    double frequence_oscillation = M_PI_2/t_end;

    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation, angle_of_attack);

    solver_all_steps(&grid);
    free_data(&grid);
    printf("M_05_moving_until_001s done");
}

// M_inf = 0.85

void M_085_static_until_001s(){
    double M_inf = 0.85;
    double t_end =  0.01;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.;
    double frequence_oscillation = 0.;

    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("M_085_static_until_001s done");

}

void M_085_moving_until_001s(){
    double M_inf = 0.85;
    double t_end =  .01;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.2; // poin
    double frequence_oscillation = M_PI_2/t_end;

    
    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("M_085_moving_until_001s done");

}


void M_085_static_until_01s(){
    double M_inf = 0.85;
    double t_end =  0.1;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.;
    double frequence_oscillation = 0.;

    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("M_085_static_until_01s done");

}

void M_085_moving_until_01s(){
    double M_inf = 0.85;
    double t_end =  .1;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.2; // poin
    double frequence_oscillation = 10.*M_PI_2/t_end;

    
    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("M_085_moving_until_01s done");

}

// M_inf > 1

void M_13_static_until_01s(){
    double M_inf = 1.3;
    double t_end =  .1;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.;
    double frequence_oscillation = 0.;

    
    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("M_13_static_until_01s done");

}

void M_13_moving_until_001s(){
    double M_inf = 1.3;
    double t_end =  .01;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.2; 
    double frequence_oscillation = 10.*M_PI_2/t_end;

    
    struct grid_and_data grid = init_grid_and_data(M_inf, grow_rate_M_inf,t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("M_13_moving_until_001s done");

}

void M_5_static_until_01s(){
    double M_inf = 5.;
    double t_end =  .1;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double grow_rate_M_inf = 0.0;

    double angle_of_attack = 0.; // en degré
    double amplitude_oscillation = 0.;
    double frequence_oscillation = 0.;

    struct grid_and_data grid = init_grid_and_data(M_inf,grow_rate_M_inf, t_end, step_max, amplitude_oscillation, frequence_oscillation,angle_of_attack);

    solver_all_steps(&grid);

    free_data(&grid);
    printf("M_5_static_until_01s done");

}

#endif 