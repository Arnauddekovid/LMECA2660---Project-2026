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

void print_values_dydi_at_j(struct grid_and_data *mygrid, int j){//ok
    printf("dydi values at j = %i \n", j);
    for(int i = 0; i<mygrid->M;i++){
        printf("%i up %f %f\n",i, mygrid->dy_updi[j][i], mygrid->y_up[j][i+1] - mygrid->y_up[j][i] );
        printf("%i down %f %f\n",i, mygrid->dy_downdi[j][i], mygrid->y_down[j][i+1] - mygrid->y_down[j][i] );
    }
}

void print_values_dydi_at_i(struct grid_and_data *mygrid, int i){
    printf("dydi values at i = %i \n", i);
    for(int j = 0; j<mygrid->N;j++){
        printf("%i up %f %f\n",j, mygrid->dy_updi[j][i], mygrid->y_up[j][i+1] - mygrid->y_up[j][i] );
        printf("%i down %f %f\n",j, mygrid->dy_downdi[j][i], mygrid->y_down[j][i+1] - mygrid->y_down[j][i] );
    }
}

void print_values_dydj_at_i(struct  grid_and_data *mygrid, int i ){
    printf("dydj values at i = %i \n", i);
    for (int j = 0; j<mygrid->N;j++){
        printf("%i up   %f %f\n", j,mygrid->dy_updj[j][i], mygrid->y_up[j+1][i]-mygrid->y_up[j][i] );
        printf("%i down %f %f\n", j,mygrid->dy_downdj[j][i], mygrid->y_down[j+1][i]-mygrid->y_down[j][i]) ;
    }
}

void print_values_dydj_at_j(struct  grid_and_data *mygrid, int j ){
    printf("dydj values at j = %i \n", j);
    for (int i = 0; i<mygrid->M;i++){
        printf("%i up   %f %f\n", i,mygrid->dy_updj[j][i], mygrid->y_up[j+1][i]-mygrid->y_up[j][i] );
        printf("%i down %f %f\n", i,mygrid->dy_downdj[j][i], mygrid->y_down[j+1][i]-mygrid->y_down[j][i]) ;
    }
}

void print_values_y(struct grid_and_data *mygrid){
    printf("y values at i = %i \n", 0);
    for (int j = 0; j < mygrid->N ; j++ ){
        printf("%i up %f down %f\n", j, mygrid->y_up[j][0], mygrid->y_down[j][0]);
    }
}

void print_x_yup_y_down_at_common_border(struct grid_and_data *mygrid){
    for (int i = 0; i<mygrid->M;i++){
        printf("index = %i , x = %f, y_up = %f, y_down = %f \n", i, mygrid->x[i], mygrid->y_up[0][i], mygrid->y_down[mygrid->N-1][i]);
    }
}

void print_leading_trailing_edges_values(struct grid_and_data *mygrid){

    printf("leading edge index = %i, trailing edge index = %i \n", mygrid->i_lead, mygrid->i_trail);
    printf("leading edge x = %f, trailing edge x = %f \n", mygrid->x[mygrid->i_lead], mygrid->x[mygrid->i_trail]);
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

void print_values_D(struct grid_and_data *mygrid, int i){
    printf("Values of D at i= %i \n",i);
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
        printf("%i up %f down %f\n", j, mygrid->D_up[j][i], mygrid->D_down[j][i]);
    }
}

void print_values_D2(struct grid_and_data *mygrid){
    for (int i = mygrid->i_lead; i <= mygrid->i_trail; i++)
    printf("i=%d D_up[0]=%g D_down[N-1]=%g dy_updj[0]=%g dy_downdj[N-1]=%g dxdi=%g\n",
        i,
        mygrid->D_up[0][i], mygrid->D_down[mygrid->N-1][i],
        mygrid->dy_updj[0][i], mygrid->dy_downdj[mygrid->N-1][i],
        mygrid->dxdi[i]);
}

void write_file_for_some_data(struct grid_and_data *mygrid,char* data_name, double** some_data_up, double** some_data_down){
    FILE * file = generate_file_other_data(mygrid,data_name);
    write_first_two_lines(mygrid, file);
    write_other_data(mygrid, file, some_data_up, some_data_down);
    fclose(file);
}

void write_data_mac_cormack(struct grid_and_data *mygrid){
    //writes two files with 0 if the j,i element has been visited bu maccormack backward and forward 
    write_file_for_some_data(mygrid, "mc f ", mygrid->mac_cormack_f_up_checker, mygrid->mac_cormack_f_down_checker);
    write_file_for_some_data(mygrid, "mc b ", mygrid->mac_cormack_b_up_checker, mygrid->mac_cormack_b_down_checker);
}

void write_data_D(struct grid_and_data *mygrid){
    //writes two files with 0 if the j,i element has been visited bu maccormack backward and forward 
    write_file_for_some_data(mygrid, "D ", mygrid->D_up, mygrid->D_down);
}

#endif