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
// gcc main.c -g -o exec -lm
// run
// bt


int main(){

    double M_inf = 0.50;
    double t_end = 0.01;
    bool airfoil = false;

    struct grid_and_data grid = init_grid_and_data(M_inf, t_end, airfoil);

    solver_all_steps(&grid);
    free_data(&grid);

}
