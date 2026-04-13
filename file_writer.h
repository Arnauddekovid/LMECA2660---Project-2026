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
    sprintf(filename,"M_inf = %f until t_end = %f.txt", infos->M_inf, infos->t_end);

    FILE * fp = fopen(filename, "w"); 
    return fp;
}

FILE * D_generate_file(struct grid_and_data *infos){
    char *filename = (char*) calloc(100, sizeof(char));
    sprintf(filename,"D M_inf = %f until t_end = %f.txt", infos->M_inf, infos->t_end);

    FILE * fp = fopen(filename, "w"); 
    return fp;
}

void write_D( struct grid_and_data *infos, FILE * fp){


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

    // write the values of D_up annd then u_down

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            fprintf(fp, "%f ", infos->D_up[infos->N-1-j][i]);
        }
        fprintf(fp, "\n");
    }
    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            fprintf(fp, "%f ", infos->D_down[infos->N-1-j][i]);

        }
        fprintf(fp, "\n");
    }
    fprintf(fp, " ");
    
}

void write_first_line( struct grid_and_data *infos, FILE *fp){

    //write x values
    for (int j= 0; j < infos->M; j++){
        fprintf(fp, "%f ", infos->x[j]);
    }
    fprintf(fp, "\n");
}

void write_one_time_step( struct grid_and_data *infos, FILE * fp){

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

    // write the values of u_up annd then u_down

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            
            double u = infos->s_up[1][infos->N-1-j][i]/infos->s_up[0][infos->N-1-j][i];
            double v = infos->s_up[2][infos->N-1-j][i]/infos->s_up[0][infos->N-1-j][i];

            double mach = sqrt(u*u+v*v)/infos->c_sound_inf;
            
            fprintf(fp, "%f ",mach);

        }
        fprintf(fp, "\n");
    }
    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            double u = infos->s_down[1][infos->N-1-j][i]/infos->s_down[0][infos->N-1-j][i];
            double v = infos->s_down[2][infos->N-1-j][i]/infos->s_down[0][infos->N-1-j][i];
            double mach = sqrt(u*u+v*v)/infos->c_sound_inf;

            fprintf(fp, "%f ", mach);
        }
        fprintf(fp, "\n");
    }

    
}


#endif