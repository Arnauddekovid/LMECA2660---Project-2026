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

void print_values_dxdi(struct grid_and_data *infos){//OK
    printf("dxdi values \n");
    for(int i = 0; i<infos->M;i++){
        printf("%i %f %f\n",i, infos->dxdi[i], infos->x[i+1] - infos->x[i]);
    }
}

void print_values_dydi_at_j(struct grid_and_data *infos, int j){//ok
    printf("dydi values at j = %i \n", j);
    for(int i = 0; i<infos->M;i++){
        printf("%i up %f %f\n",i, infos->dy_updi[j][i], infos->y_up[j][i+1] - infos->y_up[j][i] );
        printf("%i down %f %f\n",i, infos->dy_downdi[j][i], infos->y_down[j][i+1] - infos->y_down[j][i] );
    }
}

void print_values_dydi_at_i(struct grid_and_data *infos, int i){
    printf("dydi values at i = %i \n", i);
    for(int j = 0; j<infos->N;j++){
        printf("%i up %f %f\n",j, infos->dy_updi[j][i], infos->y_up[j][i+1] - infos->y_up[j][i] );
        printf("%i down %f %f\n",j, infos->dy_downdi[j][i], infos->y_down[j][i+1] - infos->y_down[j][i] );
    }
}

void print_values_dydj_at_i(struct  grid_and_data *infos, int i ){
    printf("dydj values at i = %i \n", i);
    for (int j = 0; j<infos->N;j++){
        printf("%i up   %f %f\n", j,infos->dy_updj[j][i], infos->y_up[j+1][i]-infos->y_up[j][i] );
        printf("%i down %f %f\n", j,infos->dy_downdj[j][i], infos->y_down[j+1][i]-infos->y_down[j][i]) ;
    }
}

void print_values_dydj_at_j(struct  grid_and_data *infos, int j ){
    printf("dydj values at j = %i \n", j);
    for (int i = 0; i<infos->M;i++){
        printf("%i up   %f %f\n", i,infos->dy_updj[j][i], infos->y_up[j+1][i]-infos->y_up[j][i] );
        printf("%i down %f %f\n", i,infos->dy_downdj[j][i], infos->y_down[j+1][i]-infos->y_down[j][i]) ;
    }
}

void print_values_y(struct grid_and_data *infos){
    printf("y values at i = %i \n", 0);
    for (int j = 0; j < infos->N ; j++ ){
        printf("%i up %f down %f\n", j, infos->y_up[j][0], infos->y_down[j][0]);
    }
}

void print_x_yup_y_down_at_common_border(struct grid_and_data *infos){
    for (int i = 0; i<infos->M;i++){
        printf("index = %i , x = %f, y_up = %f, y_down = %f \n", i, infos->x[i], infos->y_up[0][i], infos->y_down[infos->N-1][i]);
    }
}

void print_leading_trailing_edges_values(struct grid_and_data *infos){

    printf("leading edge index = %i, trailing edge index = %i \n", infos->i_lead, infos->i_trail);
    printf("leading edge x = %f, trailing edge x = %f \n", infos->x[infos->i_lead], infos->x[infos->i_trail]);
}

void check_conservation(struct grid_and_data *infos){
    double total_mass = 0.;
    double total_mass_star = 0.;
    for (int j = 0; j < infos->N; j++)
        for (int i = 0; i < infos->M; i++){
            total_mass += infos->s_up[0][j][i]*infos->D_up[j][i] ;
            total_mass += infos->s_down[0][j][i]*infos->D_down[j][i]  ;
            total_mass_star += infos->s_star_up[0][j][i]*infos->D_up[j][i] ;
            total_mass_star += infos->s_star_down[0][j][i]*infos->D_down[j][i];
        }
    printf("mass = %f and mass star = %f \n", total_mass, total_mass_star);
}

void print_values_D(struct grid_and_data *infos, int i){
    printf("Values of D at i= %i \n",i);
    /*
    for (int j =0; j < infos->N-1;j++){
        printf("\n values of D_up at j = %i : \n", j);
        
        for( int i = 0; i < infos->M-1; i++){
            printf("%i : %f - ", i, infos->D_up[j][i]);
        }
        printf("\n values of D_down at j = %i : \n", j);
        
        for( int i = 0; i < infos->M-1; i++){
            printf("%i : %f - ", i, infos->D_down[j][i]);
        }
    */

    for (int j =0; j <infos->N-1;j++){
        printf("%i up %f down %f\n", j, infos->D_up[j][i], infos->D_down[j][i]);
    }
}

void print_values_D2(struct grid_and_data *infos){
    for (int i = infos->i_lead; i <= infos->i_trail; i++)
    printf("i=%d D_up[0]=%g D_down[N-1]=%g dy_updj[0]=%g dy_downdj[N-1]=%g dxdi=%g\n",
        i,
        infos->D_up[0][i], infos->D_down[infos->N-1][i],
        infos->dy_updj[0][i], infos->dy_downdj[infos->N-1][i],
        infos->dxdi[i]);
}

void write_file_for_some_data(struct grid_and_data *infos,char* data_name, double** some_data_up, double** some_data_down){
    FILE * file = generate_file_other_data(infos,data_name);
    write_first_two_lines(infos, file);
    write_other_data(infos, file, some_data_up, some_data_down);
    fclose(file);
}

void write_data_D(struct grid_and_data *infos){
    //writes two files with 0 if the j,i element has been visited bu maccormack backward and forward 
    write_file_for_some_data(infos, "D ", infos->D_up, infos->D_down);
}

#endif