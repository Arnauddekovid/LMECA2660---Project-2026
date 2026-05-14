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
    int M = infos->M;
    for(int i = 0; i<infos->M;i++){
        printf("%i up %f %f\n",i, infos->dydi_up[j*M+i], infos->y_up[j*M+i+1] - infos->y_up[j*M+i] );
        printf("%i down %f %f\n",i, infos->dydi_down[j*M+i], infos->y_down[j*M+i+1] - infos->y_down[j*M+i] );
    }
}

void print_values_dydi_at_i(struct grid_and_data *infos, int i){
    printf("dydi values at i = %i \n", i);
    int M = infos->M;
    for(int j = 0; j<infos->N;j++){
        printf("%i up %f %f\n",j, infos->dydi_up[j*M+i], infos->y_up[j*M+i+1] - infos->y_up[j*M+i] );
        printf("%i down %f %f\n",j, infos->dydi_down[j*M+i], infos->y_down[j*M+i+1] - infos->y_down[j*M+i] );
    }
}

void print_values_dydj_at_i(struct  grid_and_data *infos, int i ){
    printf("dydj values at i = %i \n", i);
    int M = infos->M;
    for (int j = 0; j<infos->N;j++){
        printf("%i up   %f %f\n", j,infos->dydj_up[j*M+i], infos->y_up[(j+1)*M+i]-infos->y_up[j*M+i] );
        printf("%i down %f %f\n", j,infos->dydj_down[j*M+i], infos->y_down[(j+1)*M+i]-infos->y_down[j*M+i]) ;
    }
}

void print_values_dydj_at_j(struct  grid_and_data *infos, int j ){
    printf("dydj values at j = %i \n", j);
    int M = infos->M;
    for (int i = 0; i<infos->M;i++){
        printf("%i up   %f %f\n", i,infos->dydj_up[j*M+i], infos->y_up[(j+1)*M+i]-infos->y_up[j*M+i] );
        printf("%i down %f %f\n", i,infos->dydj_down[j*M+i], infos->y_down[(j+1)*M+i]-infos->y_down[j*M+i]) ;
    }
}

void print_values_y(struct grid_and_data *infos){
    printf("y values at i = %i \n", 0);
    int M = infos->M;
    for (int j = 0; j < infos->N ; j++ ){
        printf("%i up %f down %f\n", j, infos->y_up[j*M], infos->y_down[j*M]);
    }
}

void print_x_yup_y_down_at_common_border(struct grid_and_data *infos){
    int M = infos->M;

    for (int i = 0; i<infos->M;i++){
        printf("index = %i , x = %f, y_up = %f, y_down = %f \n", i, infos->x[i], infos->y_up[i], infos->y_down[(infos->N-1)*M+i]);
    }
}

void print_leading_trailing_edges_values(struct grid_and_data *infos){

    printf("leading edge index = %i, trailing edge index = %i \n", infos->i_lead, infos->i_trail);
    printf("leading edge x = %f, trailing edge x = %f \n", infos->x[infos->i_lead], infos->x[infos->i_trail]);
}

void check_conservation(struct grid_and_data *infos){
    double total_mass = 0.;
    double total_mass_star = 0.;
    int M = infos->M;
    for (int j = 0; j < infos->N; j++)
        for (int i = 0; i < infos->M; i++){
            total_mass += infos->s_up[j*M+i]*infos->D_up[j*M+i] ;
            total_mass += infos->s_down[j*M+i]*infos->D_down[j*M+i]  ;
            total_mass_star += infos->s_star_up[j*M+i]*infos->D_up[j*M+i] ;
            total_mass_star += infos->s_star_down[j*M+i]*infos->D_down[j*M+i];
        }
    printf("mass = %f and mass star = %f \n", total_mass, total_mass_star);
}

void print_values_D(struct grid_and_data *infos, int i){
    printf("Values of D at i= %i \n",i);
    int M = infos->M;

    for (int j =0; j <infos->N-1;j++){
        printf("%i up %f down %f\n", j, infos->D_up[j*M+i], infos->D_down[j*M+i]);
    }
}

void print_values_D2(struct grid_and_data *infos){
    int M = infos->M;
    for (int i = infos->i_lead; i <= infos->i_trail; i++)
    printf("i=%d D_up[0]=%g D_down[N-1]=%g dydj_up[0]=%g dydj_down[N-1]=%g dxdi=%g\n",
        i,
        infos->D_up[i], infos->D_down[(infos->N-1)*M+i],
        infos->dydj_up[i], infos->dydj_down[(infos->N-1)*M+i],
        infos->dxdi[i]);
}






#endif