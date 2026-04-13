//
//  solver.h
//  MECA2660 projet 2026


#ifndef debug_functions_h
#define debug_functions_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  
#include "grid_and_data.h"
#include "file_writer.h"

void print_values_dxdi(struct grid_and_data *mygrid){//OK
    printf("dxdi values \n");
    for(int i = 0; i<mygrid->M;i++){
        printf("%i %f %f\n",i, mygrid->dxdi[i], mygrid->x[i+1] - mygrid->x[i]);
    }
}

void print_values_dydi(struct grid_and_data *mygrid, int j){//ok
    printf("dydi values at j = %i \n", j);
    for(int i = 0; i<mygrid->M;i++){
        printf("%i up %f %f\n",i, mygrid->dy_updi[j][i], mygrid->y_up[j][i+1] - mygrid->y_up[j][i] );
        printf("%i down %f %f\n",i, mygrid->dy_downdi[j][i], mygrid->y_down[j][i+1] - mygrid->y_down[j][i] );
    }
    
}

void print_values_dydj(struct  grid_and_data *mygrid, int i ){
    printf("dydj values at i = %i \n", i);
    for (int j = 0; j<mygrid->N;j++){
        printf("%i up   %f %f\n", j,mygrid->dy_updj[j][i], mygrid->y_up[j+1][i]-mygrid->y_up[j][i] );
        printf("%i down %f %f\n", j,mygrid->dy_downdj[j][i], mygrid->y_down[j+1][i]-mygrid->y_down[j][i]) ;
    }
}

void print_values_y(struct grid_and_data *mygrid){
    printf("y values at i = %i \n", 0);
    for (int j = 0; j < mygrid->N ; j++ ){
        printf("%i up %f down %f\n", j, mygrid->y_up[j][0], mygrid->y_down[j][0]);
    }
}

void print_values_D(struct grid_and_data *mygrid){
    printf("Values of D \n");
    /*
    for (int j =0; j < mygrid->N-1;j++){
        printf("\n values of D_up at j = %i : \n", j);
        
        for( int i = 0; i < mygrid->M-1; i++){
            printf("%i : %f - ", i, mygrid->D_up[j][i]);
        }
        printf("\n values of D_down at j = %i : \n", j);
        
        for( int i = 0; i < mygrid->M-1; i++){
            printf("%i : %f - ", i, mygrid->D_down[j][i]);
        }
    */

    for (int j =0; j <mygrid->N-1;j++){
        printf("%i up %f down %f\n", j, mygrid->D_up[j][0], mygrid->D_down[j][0]);
    }
}

void write_D_file(struct grid_and_data *mygrid){
    FILE * file = D_generate_file(mygrid);
    write_first_line(mygrid, file);
    write_D(mygrid, file);
    fclose(file);
}

#endif