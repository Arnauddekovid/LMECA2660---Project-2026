#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  
#include <complex.h>

#include "grid_and_data.h"
#include "file_writer.h"
#include "solver.h"

#include "debug_functions.h"

// move de bac1 je sais
// cd Projet
// gcc main.c -o exec
// ./exec

// trouver les segfaults 
// gcc main.c -g
// gcc main.c -g -o exec -lmx
// run
// btgcc main.c -o exec

int main(){

    double M_inf = 0.85;
    double t_end = .1;
    double dt = 1e-5;

    double step_max_double = t_end*1000.;// tt les combien de step on stocke les valeurs dans le fichier ?, ave ct_end*1000 et dt = 1e-5 on a une animation de 3 sec environ
    int step_max = (int) step_max_double ;

    double C = 0.01; //on met à 0 pr M_inf = 0.5
    double eps_min = 2.; // pour M_inf = 0.5 , on prend eps_min = 0.1

    struct grid_and_data grid = init_grid_and_data(M_inf, t_end, dt, step_max, C, eps_min);

    solver_all_steps(&grid);

    //write_data_D(&grid);

    free_data(&grid);
    printf("done with M_inf = %f, eps_min = %f, C = %f ",M_inf, eps_min, C);

    return 0;

}
