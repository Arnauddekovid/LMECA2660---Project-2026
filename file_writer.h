//
//  file_writer.h
//  MECA2660 projet 2026

#ifndef file_writer_h
#define file_writer_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  

#include "grid_and_data.h"

FILE * generate_file(struct grid_and_data *infos){
    char *filename = (char*) calloc(100, sizeof(char));
    if (infos->mu > 0.){
        sprintf(filename,"M_%g_moving_until_%gs.txt", infos->M_inf,infos->t_end);
    }
    else if (infos->grow_rate_M_inf > 0.){
        sprintf(filename,"M_%g_growing_until_%gs.txt", infos->M_inf,infos->t_end);  
    }
    else{
        sprintf(filename,"M_%g_static_until_%gs.txt", infos->M_inf,infos->t_end);
    }
    FILE * fp = fopen(filename, "w"); 
    return fp;
}

FILE * generate_file_other_data(struct grid_and_data *infos, char * data_name){
    char *filename = (char*) calloc(100, sizeof(char));
    sprintf(filename,"%s M_inf = %f until t_end = %f.txt",data_name, infos->M_inf, infos->t_end);
    FILE * fp = fopen(filename, "w"); 
    return fp;
}

void write_other_data( struct grid_and_data *infos, FILE * fp, double** data_up, double** data_down){

    // write the value of t
    fprintf(fp, "%f \n", infos->t);
    //write the coordinates of y_up
    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            fprintf(fp, "%f ", infos->y_up[infos->N-1-j][i]);
        }
        fprintf(fp, "\n");
    }
    //write the coordinates of y_down

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            fprintf(fp, "%f ", infos->y_down[infos->N-1-j][i]);
        }
        fprintf(fp, "\n");
    }
    // rewrite the value of t bc why not
    fprintf(fp, "%f \n", infos->t);

    // write the values of data up annd then data down

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            fprintf(fp, "%f ", data_up[infos->N-1-j][i]);
        }
        fprintf(fp, "\n");
    }
    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            fprintf(fp, "%f ", data_down[infos->N-1-j][i]);

        }
        fprintf(fp, "\n");
    }
    fprintf(fp, " ");
    
}

void write_first_two_lines( struct grid_and_data *infos, FILE *fp){

    fprintf(fp, "M_inf = %g T_0_inf = %g p_0_inf = %g R_fluid = %g gamma = %g c = %g mu = %g omega = %g \n", infos->M_inf, infos->T_0_inf, infos->p_0_inf, infos->R_fluid, infos->gamma, infos->c, infos->mu, infos->omega);
    for (int j= 0; j < infos->M; j++){
        fprintf(fp, "%f ", infos->x[j]);
    }
    fprintf(fp, "\n");
}



void write_one_time_step( struct grid_and_data *infos, FILE * fp){

    // write the value of t
    fprintf(fp, "y_grid at t = %g M_inf = %g\n", infos->t, infos->M_inf  );
    //write the coordinates of y_up
    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            fprintf(fp, "%f ", infos->y_up[infos->N-1-j][i]);
        }
        fprintf(fp, "\n");
    }
    //write the coordinates of y_down

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            fprintf(fp, "%f ", infos->y_down[infos->N-1-j][i]);
        }
        fprintf(fp, "\n");
    }

    for (int k = 0; k < 4; k++){
        // rewrite the value of t bc why not
        fprintf(fp, "s_%i at t = %g M_inf = %g\n",k, infos->t, infos->M_inf   );
        // write the values of s_0 up annd then s_0 down
        for (int j = 0; j < infos->N; j++){
            for (int i = 0; i < infos->M; i++){
                fprintf(fp, "%f ", infos->s_up[k][infos->N-1-j][i]);
            }
            fprintf(fp, "\n");
        }
        for (int j = 0; j < infos->N; j++){
            for (int i = 0; i < infos->M; i++){
                fprintf(fp, "%f ",  infos->s_down[k][infos->N-1-j][i]);
            }
            fprintf(fp, "\n");
        }
    }
}

#endif