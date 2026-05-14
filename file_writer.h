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


char * generate_filename(struct grid_and_data *infos){
    /*
        pre :
            infos : the structure we use
        post :
           char * : generates binary file name with Resolution_M_inf__moving_or_not_until_t_end_binary.bin and store in the hdd or not
    */
    char *filename = (char*) calloc(170, sizeof(char));
    char *resolution = (char*) calloc(10, sizeof(char));
    if (infos->resolution == 0.){
        snprintf(resolution,10,"low_res_");
    }
    else if (infos->resolution == 0.5){
        snprintf(resolution,10,"mid_res_");
    }
    else if (infos->resolution == 1.){
        snprintf(resolution,10,"high_res_");
    }
    else{
        snprintf(resolution,10,"low_res_" );

    }

    char *M_inf = (char*) calloc(20, sizeof(char));
    snprintf(M_inf,20,"M_inf_%g_", infos->M_inf);

    char* behavior = (char*) calloc(25, sizeof(char));
    if (infos->grow_rate_M_inf != 0.){
        snprintf(behavior,25,"growing_rate_M_inf_%g_", infos->grow_rate_M_inf);
    }
    else if (infos->mu != 0.){
        snprintf(behavior,25,"oscillating_mu_%g_", infos->mu);  
    }
    else {
        double AOA_deg = infos->AOA*180./M_PI;
        snprintf(behavior,25,"static_A0A_%g_", AOA_deg);
    }

    char *t_end = (char*) calloc(20, sizeof(char));
    snprintf(t_end,20,"until_%gs", infos->t_end);

    if (infos->hdd_storage){
        snprintf(filename,170,"/Volumes/Bruno/binaryFiles/%s%s%s%s.bin", resolution, M_inf, behavior, t_end);
    }
    else{
        snprintf(filename,170,"binaryFiles/%s%s%s%s.bin", resolution, M_inf, behavior, t_end);
    }

    return filename;
}

void write_first_two_lines(struct grid_and_data *infos, FILE *fp) {
    /*
        pre :
            infos : the structure we use
            fp : the binary file where
        post :
           writes the header of the file with some informations and then the x 1D array
    */

    // Write header as raw binary
    fwrite(&infos->N, sizeof(int), 1, fp);
    fwrite(&infos->M, sizeof(int), 1, fp);
    fwrite(&infos->M_inf, sizeof(double), 1, fp);
    fwrite(&infos->T_0_inf, sizeof(double), 1, fp);
    fwrite(&infos->p_0_inf, sizeof(double), 1, fp);
    fwrite(&infos->R_fluid, sizeof(double), 1, fp);
    fwrite(&infos->gamma, sizeof(double), 1, fp);
    fwrite(&infos->c, sizeof(double), 1, fp);
    fwrite(&infos->mu, sizeof(double), 1, fp);
    fwrite(&infos->omega, sizeof(double), 1, fp);
    // Write x grid (M values)
    fwrite(infos->x, sizeof(double), infos->M, fp);
}

void write_one_time_step(struct grid_and_data *infos, FILE *fp) {
    /*
        pre :
            infos : the structure we use
            fp : the binary file where
        post :
           writes the timestep, the Mach Inf, the y grid and the characteristics 
    */
    int M = infos->M;
    int N = infos->N;
    
    fwrite(&infos->t, sizeof(double), 1, fp);
    fwrite(&infos->M_inf, sizeof(double), 1, fp);

    fwrite(infos->y_up, sizeof(double), N*M, fp);
    fwrite(infos->y_down, sizeof(double), N*M, fp);

    fwrite(infos->s_up, sizeof(double), 4*N*M, fp);
    fwrite(infos->s_down, sizeof(double), 4*N*M, fp);
    
}


#endif