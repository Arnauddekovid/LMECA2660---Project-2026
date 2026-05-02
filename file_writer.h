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

    //write x valuesx
    fprintf(fp, "eps_min = %f, C = %f, \n", infos->eps_min, infos->C);
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

    double gamma = infos->gamma;
    double T_0 = infos->T_0_inf; 
    double R = infos->R_fluid;

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            double rho = infos->s_up[0][infos->N-1-j][i];
            double u = infos->s_up[1][infos->N-1-j][i]/rho;
            double v = infos->s_up[2][infos->N-1-j][i]/rho;
            double U0 = infos->s_up[3][infos->N-1-j][i]/rho;
            double p = rho*(gamma-1.)*(U0 -0.5*(u*u+v*v));

            double c_sound  = sqrt(gamma*p/rho);

            double T = p/(rho*R);

            //double mach = sqrt(2./(gamma-1.)*(T_0/T-1.));
            double mach = sqrt(u*u+v*v)/c_sound;
            double value =  mach;
            
            fprintf(fp, "%f ", value);

        }
        fprintf(fp, "\n");
    }
    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            double rho = infos->s_down[0][infos->N-1-j][i];
            double u = infos->s_down[1][infos->N-1-j][i]/rho;
            double v = infos->s_down[2][infos->N-1-j][i]/rho;
            double U0 =infos->s_down[3][infos->N-1-j][i]/rho;
            double p = rho*(gamma-1.)*(U0 -0.5*(u*u+v*v));

            double c_sound  = sqrt(gamma*p/rho);

            double T = p/(rho*R);

            //double mach = sqrt(2./(gamma-1.)*(T_0/T-1.));
            double mach = sqrt(u*u+v*v)/c_sound;
            double value = mach;
            
            fprintf(fp, "%f ", value);
        }
        fprintf(fp, "\n");
    }
    
}

#endif