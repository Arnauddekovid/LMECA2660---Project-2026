//
//  solver.h
//  MECA2660 projet 2026

#ifndef solver_h
#define solver_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  
#include <sys/stat.h>
#include <sys/types.h>

#include "grid_and_data.h"



typedef void (* bounds_s_function)(struct grid_and_data *);

void set_init_s(struct grid_and_data *infos){
    int N = infos->N;
    int M = infos->M;
    double gamma  = infos->gamma;
    double rho_inf = infos->rho_inf;
    double p_inf   = infos->p_inf;
    double u_inf   = infos->M_inf * infos->c_sound_inf;
    double v_inf = 0.;
    double U0_inf = p_inf / ((gamma - 1.) * rho_inf) + 0.5 * u_inf * u_inf;      

    for (int j = 0; j < N; j++){
        for (int i = 0; i < M; i++){
            infos->s_up[0*M*N+j*M+i] = rho_inf;
            infos->s_up[1*M*N+j*M+i] = rho_inf * u_inf;
            infos->s_up[2*M*N+j*M+i] = 0.;
            infos->s_up[3*M*N+j*M+i] = rho_inf * U0_inf;

            infos->s_down[0*M*N+j*M+i] = rho_inf;
            infos->s_down[1*M*N+j*M+i] = rho_inf * u_inf;
            infos->s_down[2*M*N+j*M+i] = 0.;
            infos->s_down[3*M*N+j*M+i] = rho_inf * U0_inf;

            for (int k = 0; k < 4; k++){
                infos->s_star_up[k*M*N+j*M+i]   = infos->s_up[k*M*N+j*M+i];
                infos->s_star_down[k*M*N+j*M+i] = infos->s_down[k*M*N+j*M+i];
            }
        }
    }

    
}

// subsonic case

void set_bounds_s_subsonic(struct grid_and_data *infos ){

    int N = infos->N;
    int M = infos->M;
    double rho_inf = infos->rho_inf;
    double p_inf = infos->p_inf;
    double u_inf = infos->M_inf*infos->c_sound_inf;
    double v_inf = 0.;
    double v_0 = 0.;

    double rho_0;
    double p_0;
    double u_0;
    double arg;
    double U_0;

    double rho_end;
    double u_end;
    double v_end;
    double U0_end;


    for (int j = 0; j< N; j++){
        //inlet conditions

        rho_0 = infos->s_up[j*M+1]; 
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        arg = u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.) * (p_inf/rho_inf - p_0/rho_0);
        u_0 = sqrt(fabs(arg));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);

        infos->s_up[j*M] = rho_0;
        infos->s_up[1*N*M+j*M] = rho_0*u_0;
        infos->s_up[2*N*M+j*M] = rho_0*v_0;
        infos->s_up[3*N*M+j*M] = rho_0*U_0;

        rho_0 = infos->s_down[j*M+1];
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        arg = u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.) * (p_inf/rho_inf - p_0/rho_0);
        u_0 = sqrt(fabs(arg));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);

        infos->s_down[j*M]= rho_0;
        infos->s_down[1*N*M+j*M]= rho_0*u_0;
        infos->s_down[2*N*M+j*M]= rho_0*v_0;
        infos->s_down[3*N*M+j*M]= rho_0*U_0;

        rho_end = infos->s_up[0*N*M+j*M+M-2];
        u_end = infos->s_up[1*N*M+j*M+M-2]/rho_end;
        v_end = infos->s_up[2*N*M+j*M+M-2]/rho_end;

        U0_end = p_inf/((infos->gamma-1.)*rho_end) + 0.5*(u_end*u_end+v_end*v_end);

        infos->s_up[0*N*M+j*M+M-1] = infos->s_up[0*N*M+j*M+M-2];
        infos->s_up[1*N*M+j*M+M-1] = infos->s_up[1*N*M+j*M+M-2];
        infos->s_up[2*N*M+j*M+M-1] = infos->s_up[2*N*M+j*M+M-2];
        infos->s_up[3*N*M+j*M+M-1] = rho_end*U0_end;

        rho_end = infos->s_down[0*N*M+j*M+M-2];
        u_end = infos->s_down[1*N*M+j*M+M-2]/rho_end;
        v_end = infos->s_down[2*N*M+j*M+M-2]/rho_end;
        U0_end = p_inf/((infos->gamma-1.)*rho_end) + 0.5*(u_end*u_end+v_end*v_end);

        infos->s_down[0*N*M+j*M+M-1] = infos->s_down[0*N*M+j*M+M-2];
        infos->s_down[1*N*M+j*M+M-1] = infos->s_down[1*N*M+j*M+M-2];
        infos->s_down[2*N*M+j*M+M-1] = infos->s_down[2*N*M+j*M+M-2];
        infos->s_down[3*N*M+j*M+M-1] = rho_end*U0_end;
    }

    // wall dumby bumd conditions ( on prend pas encore en compte les trucs courbés)
    // à priori c'es u dx - v dy = 0 mais à vérifier

    double rho_wall; 
    double u_wall;
    double v_wall;
    double U_wall;

    double c_y0;
    double dx;
    double dy;
    double norm_direction;

    // wall

    for (int i = 1; i <M-1; i++){

        infos->s_up[0*N*M+(N-1)*M+i] = infos->s_up[0*N*M+(N-2)*M+i];
        infos->s_up[1*N*M+(N-1)*M+i] = infos->s_up[1*N*M+(N-2)*M+i];
        infos->s_up[2*N*M+(N-1)*M+i] = 0.0;
        infos->s_up[3*N*M+(N-1)*M+i] = infos->s_up[3*N*M+(N-2)*M+i];

        infos->s_down[0*N*M+0*M+i] = infos->s_down[0*N*M+1*M+i];
        infos->s_down[1*N*M+0*M+i] = infos->s_down[1*N*M+1*M+i];
        infos->s_down[2*N*M+0*M+i] = 0.0;
        infos->s_down[3*N*M+0*M+i] = infos->s_down[3*N*M+1*M+i];
    }

    //upstream
    
    for (int i = 0; i<infos->i_lead;i++){
        for (int k = 0; k<4;k++){
            infos->s_up[k*N*M+0*M+i] = infos->s_down[k*N*M+(N-1)*M+i];
        }
    }

    //airfoil
    for(int i = infos->i_lead; i < infos->i_trail+1; i++){
        // up : j=0
        double rho = infos->s_up[0*N*M+1*M+i];
        double u   = infos->s_up[1*N*M+1*M+i] / rho;
        double v   = infos->s_up[2*N*M+1*M+i] / rho;
        
        // normale à la paroi (perpendiculaire à la tangente (dx, dy_updi))
        double dx = infos->dxdi[i];
        double dy = infos->dydi_up[0*M+i];
        double norm2 = dx*dx + dy*dy;

        double c_y0 = (u*dx + v*dy)/norm2;

        infos->s_up[0*N*M+0*M+i] = rho;
        infos->s_up[1*N*M+0*M+i] = rho *c_y0*dx ;
        infos->s_up[2*N*M+0*M+i] = rho *c_y0*dy ;
        infos->s_up[3*N*M+0*M+i] = infos->s_up[3*N*M+1*M+i];

        // down : j=N-1, voisin = j=N-2
        rho = infos->s_down[0*N*M+(N-2)*M+i];
        u   = infos->s_down[1*N*M+(N-2)*M+i] / rho;
        v   = infos->s_down[2*N*M+(N-2)*M+i] / rho;
        
        dy = infos->dydi_down[(N-1)*M+i];
        norm2 = dx*dx + dy*dy;
        c_y0 = (u*dx + v*dy)/norm2;
        
        infos->s_down[0*N*M+(N-1)*M+i] = rho;
        infos->s_down[1*N*M+(N-1)*M+i] = rho *c_y0*dx ;
        infos->s_down[2*N*M+(N-1)*M+i] = rho *c_y0*dy ;
        infos->s_down[3*N*M+(N-1)*M+i] = infos->s_down[3*N*M+(N-2)*M+i];
    }


    for (int i = infos->i_trail+1; i<M;i++){
        for (int k = 0; k<4;k++){
            infos->s_up[k*N*M+0*M+i] = infos->s_down[k*N*M+(N-1)*M+i];
        }
    }

}

void set_bounds_s_star_subsonic(struct grid_and_data *infos ){
    int N = infos->N;
    int M = infos->M;
    double rho_inf = infos->rho_inf;
    double p_inf = infos->p_inf;
    double u_inf = infos->M_inf*infos->c_sound_inf;
    double v_0 = 0;

    double rho_0, p_0, u_0, U_0;
    double rho_end, u_end, v_end, U0_end;


    for (int j = 0; j < N; j++){
        // inlet conditions
        rho_0 = infos->s_star_up[j*M+1];
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        u_0 = sqrt(u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.)*(p_inf/rho_inf - p_0/rho_0));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);
        infos->s_star_up[j*M] = rho_0;
        infos->s_star_up[1*N*M+j*M] = rho_0*u_0;
        infos->s_star_up[2*N*M+j*M] = rho_0*v_0;
        infos->s_star_up[3*N*M+j*M] = rho_0*U_0;

        rho_0 = infos->s_star_down[j*M+1];
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        u_0 = sqrt(u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.)*(p_inf/rho_inf - p_0/rho_0));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);
        infos->s_star_down[j*M] = rho_0;
        infos->s_star_down[1*N*M+j*M] = rho_0*u_0;
        infos->s_star_down[2*N*M+j*M] = rho_0*v_0;
        infos->s_star_down[3*N*M+j*M] = rho_0*U_0;
    }

    for (int j = 0; j < N; j++){
        // outlet conditions
        rho_end = infos->s_star_up[0*N*M+j*M+M-2];
        u_end = infos->s_star_up[1*N*M+j*M+M-2]/rho_end;
        v_end = infos->s_star_up[2*N*M+j*M+M-2]/rho_end;
        U0_end = p_inf/((infos->gamma-1.)*rho_end) + 0.5*(u_end*u_end+v_end*v_end);
        infos->s_star_up[0*N*M+j*M+M-1] = infos->s_star_up[0*N*M+j*M+M-2];
        infos->s_star_up[1*N*M+j*M+M-1] = infos->s_star_up[1*N*M+j*M+M-2];
        infos->s_star_up[2*N*M+j*M+M-1] = infos->s_star_up[2*N*M+j*M+M-2];
        infos->s_star_up[3*N*M+j*M+M-1] = rho_end*U0_end;

        rho_end = infos->s_star_down[0*N*M+j*M+M-2];
        u_end = infos->s_star_down[1*N*M+j*M+M-2]/rho_end;
        v_end = infos->s_star_down[2*N*M+j*M+M-2]/rho_end;
        U0_end = p_inf/((infos->gamma-1.)*rho_end) + 0.5*(u_end*u_end+v_end*v_end);
        infos->s_star_down[0*N*M+j*M+M-1] = infos->s_star_down[0*N*M+j*M+M-2];
        infos->s_star_down[1*N*M+j*M+M-1] = infos->s_star_down[1*N*M+j*M+M-2];
        infos->s_star_down[2*N*M+j*M+M-1] = infos->s_star_down[2*N*M+j*M+M-2];
        infos->s_star_down[3*N*M+j*M+M-1] = rho_end*U0_end;
    }

    //  walls
    for (int i = 1; i < M-1; i++){
        infos->s_star_up[0*N*M+(N-1)*M+i] = infos->s_star_up[0*N*M+(N-2)*M+i];
        infos->s_star_up[1*N*M+(N-1)*M+i] = infos->s_star_up[1*N*M+(N-2)*M+i];
        infos->s_star_up[2*N*M+(N-1)*M+i] = 0.;
        infos->s_star_up[3*N*M+(N-1)*M+i] = infos->s_star_up[3*N*M+(N-2)*M+i];

        infos->s_star_down[0*N*M+0*M+i] = infos->s_star_down[0*N*M+1*M+i];
        infos->s_star_down[1*N*M+0*M+i] = infos->s_star_down[1*N*M+1*M+i];
        infos->s_star_down[2*N*M+0*M+i] = 0.;
        infos->s_star_down[3*N*M+0*M+i] = infos->s_star_down[3*N*M+1*M+i];
    }

    // frontière commune upstream
    for (int i = 0; i < infos->i_lead; i++){
        for (int k = 0; k < 4; k++)
            infos->s_star_up[k*N*M+0*M+i] = infos->s_star_down[k*N*M+(N-1)*M+i];
    }

    //airfoil
    for(int i = infos->i_lead; i < infos->i_trail; i++){
        // up : j=0
        double rho = infos->s_star_up[0*N*M+1*M+i];
        double u   = infos->s_star_up[1*N*M+1*M+i] / rho;
        double v   = infos->s_star_up[2*N*M+1*M+i] / rho;
        
        // normale à la paroi (perpendiculaire à la tangente (dx, dy_updi))
        double dx = infos->dxdi[i];
        double dy = infos->dydi_up[0*M+i];
        double norm2 = dx*dx + dy*dy;

        double c_y0 = (u*dx + v*dy)/norm2;

        infos->s_star_up[0*N*M+0*M+i] = rho;
        infos->s_star_up[1*N*M+0*M+i] = rho *c_y0*dx ;
        infos->s_star_up[2*N*M+0*M+i] = rho *c_y0*dy ;
        infos->s_star_up[3*N*M+0*M+i] = infos->s_star_up[3*N*M+1*M+i];

        // down : j=N-1, voisin = j=N-2
        rho = infos->s_star_down[0*N*M+(N-2)*M+i];
        u   = infos->s_star_down[1*N*M+(N-2)*M+i] / rho;
        v   = infos->s_star_down[2*N*M+(N-2)*M+i] / rho;
        
        dy = infos->dydi_down[(N-1)*M+i];
        norm2 = dx*dx + dy*dy;
        c_y0 = (u*dx + v*dy)/norm2;
        
        infos->s_star_down[0*N*M+(N-1)*M+i] = rho;
        infos->s_star_down[1*N*M+(N-1)*M+i] = rho *c_y0*dx ;
        infos->s_star_down[2*N*M+(N-1)*M+i] = rho *c_y0*dy ;
        infos->s_star_down[3*N*M+(N-1)*M+i] = infos->s_star_down[3*N*M+(N-2)*M+i];
    }

    // frontière commune downstream
    for (int i = infos->i_trail+1; i < M; i++){
        for (int k = 0; k < 4; k++)
            infos->s_star_up[k*N*M+0*M+i] = infos->s_star_down[k*N*M+(N-1)*M+i];
    }
}

// supersonic case

void set_bounds_s_supersonic(struct grid_and_data *infos ){
    int N = infos->N;
    int M = infos->M;
    double rho_inf = infos->rho_inf;
    double p_inf = infos->p_inf;
    double u_inf = infos->M_inf*infos->c_sound_inf;
    double v_inf = 0.;
    double v_0 = 0.;

    double rho_0;
    double p_0;
    double u_0;
    double arg;
    double U_0;

    double rho_end;
    double u_end;
    double v_end;
    double U0_end;


    for (int j = 0; j< N; j++){
        //inlet

        U_0 = p_inf/((infos->gamma-1.)*rho_inf) + 0.5*(u_inf*u_inf);
        infos->s_up[j*M] = rho_inf;
        infos->s_up[1*N*M+j*M] = rho_inf*u_inf;
        infos->s_up[2*N*M+j*M] = 0.;
        infos->s_up[3*N*M+j*M] = rho_inf*U_0;


        infos->s_down[j*M]= rho_inf;
        infos->s_down[1*N*M+j*M]= rho_inf*u_inf;
        infos->s_down[2*N*M+j*M]= 0.;
        infos->s_down[3*N*M+j*M]= rho_inf*U_0;


        //outlet condintions

        infos->s_up[0*N*M+j*M+M-1] = infos->s_up[0*N*M+j*M+M-2];
        infos->s_up[1*N*M+j*M+M-1] = infos->s_up[1*N*M+j*M+M-2];
        infos->s_up[2*N*M+j*M+M-1] = infos->s_up[2*N*M+j*M+M-2];
        infos->s_up[3*N*M+j*M+M-1] = infos->s_up[3*N*M+j*M+M-2];


        infos->s_down[0*N*M+j*M+M-1] = infos->s_down[0*N*M+j*M+M-2];
        infos->s_down[1*N*M+j*M+M-1] = infos->s_down[1*N*M+j*M+M-2];
        infos->s_down[2*N*M+j*M+M-1] = infos->s_down[2*N*M+j*M+M-2];
        infos->s_down[3*N*M+j*M+M-1] = infos->s_down[3*N*M+j*M+M-2];
    }


    // wall dumby bumd conditions ( on prend pas encore en compte les trucs courbés)
    // à priori c'es u dx - v dy = 0 mais à vérifier

    double rho_wall; 
    double u_wall;
    double v_wall;
    double U_wall;

    double c_y0;
    double dx;
    double dy;
    double norm_direction;

    // wall

    for (int i = 1; i <M-1; i++){

        infos->s_up[0*N*M+(N-1)*M+i] = infos->s_up[0*N*M+(N-2)*M+i];
        infos->s_up[1*N*M+(N-1)*M+i] = infos->s_up[1*N*M+(N-2)*M+i];
        infos->s_up[2*N*M+(N-1)*M+i] = 0.0;
        infos->s_up[3*N*M+(N-1)*M+i] = infos->s_up[3*N*M+(N-2)*M+i];

        infos->s_down[0*N*M+0*M+i] = infos->s_down[0*N*M+1*M+i];
        infos->s_down[1*N*M+0*M+i] = infos->s_down[1*N*M+1*M+i];
        infos->s_down[2*N*M+0*M+i] = 0.0;
        infos->s_down[3*N*M+0*M+i] = infos->s_down[3*N*M+1*M+i];
    }

    //upstream
    
    for (int i = 0; i<infos->i_lead;i++){
        for (int k = 0; k<4;k++){
            infos->s_up[k*N*M+0*M+i] = infos->s_down[k*N*M+(N-1)*M+i];
        }
    }

    //airfoil
    for(int i = infos->i_lead; i < infos->i_trail+1; i++){
        // up : j=0
        double rho = infos->s_up[0*N*M+1*M+i];
        double u   = infos->s_up[1*N*M+1*M+i] / rho;
        double v   = infos->s_up[2*N*M+1*M+i] / rho;
        
        // normale à la paroi (perpendiculaire à la tangente (dx, dy_updi))
        double dx = infos->dxdi[i];
        double dy = infos->dydi_up[0*M+i];
        double norm2 = dx*dx + dy*dy;

        double c_y0 = (u*dx + v*dy)/norm2;

        infos->s_up[0*N*M+0*M+i] = rho;
        infos->s_up[1*N*M+0*M+i] = rho *c_y0*dx ;
        infos->s_up[2*N*M+0*M+i] = rho *c_y0*dy ;
        infos->s_up[3*N*M+0*M+i] = infos->s_up[3*N*M+1*M+i];

        // down : j=N-1, voisin = j=N-2
        rho = infos->s_down[0*N*M+(N-2)*M+i];
        u   = infos->s_down[1*N*M+(N-2)*M+i] / rho;
        v   = infos->s_down[2*N*M+(N-2)*M+i] / rho;
        
        dy = infos->dydi_down[(N-1)*M+i];
        norm2 = dx*dx + dy*dy;
        c_y0 = (u*dx + v*dy)/norm2;
        
        infos->s_down[0*N*M+(N-1)*M+i] = rho;
        infos->s_down[1*N*M+(N-1)*M+i] = rho *c_y0*dx ;
        infos->s_down[2*N*M+(N-1)*M+i] = rho *c_y0*dy ;
        infos->s_down[3*N*M+(N-1)*M+i] = infos->s_down[3*N*M+(N-2)*M+i];
    }


    for (int i = infos->i_trail+1; i<M;i++){
        for (int k = 0; k<4;k++){
            infos->s_up[k*N*M+0*M+i] = infos->s_down[k*N*M+(N-1)*M+i];
        }
    }

}

void set_bounds_s_star_supersonic(struct grid_and_data *infos ){
    int N = infos->N;
    int M = infos->M;
    double rho_inf = infos->rho_inf;
    double p_inf = infos->p_inf;
    
    double u_inf = infos->M_inf*infos->c_sound_inf;
    double v_0 = 0;

    double rho_0, p_0, u_0, U_0;
    double rho_end, u_end, v_end, U0_end;


    for (int j = 0; j < N; j++){
        // inlet conditions on met tout en inf 
        rho_0 = rho_inf;
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        u_0 = sqrt(u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.)*(p_inf/rho_inf - p_0/rho_0));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);
        infos->s_star_up[j*M] = rho_0;
        infos->s_star_up[1*N*M+j*M] = rho_0*u_0;
        infos->s_star_up[2*N*M+j*M] = rho_0*v_0;
        infos->s_star_up[3*N*M+j*M] = rho_0*U_0;

        rho_0 = rho_inf;
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        u_0 = sqrt(u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.)*(p_inf/rho_inf - p_0/rho_0));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);
        infos->s_star_down[j*M] = rho_0;
        infos->s_star_down[1*N*M+j*M] = rho_0*u_0;
        infos->s_star_down[2*N*M+j*M] = rho_0*v_0;
        infos->s_star_down[3*N*M+j*M] = rho_0*U_0;
    }

    for (int j = 0; j < N; j++){
        // outlet conditions on extrapole 

        infos->s_star_up[0*N*M+j*M+M-1] = infos->s_star_up[0*N*M+j*M+M-2];
        infos->s_star_up[1*N*M+j*M+M-1] = infos->s_star_up[1*N*M+j*M+M-2];
        infos->s_star_up[2*N*M+j*M+M-1] = infos->s_star_up[2*N*M+j*M+M-2];
        infos->s_star_up[3*N*M+j*M+M-1] = infos->s_star_up[3*N*M+j*M+M-2];

        infos->s_star_down[0*N*M+j*M+M-1] = infos->s_star_down[0*N*M+j*M+M-2];
        infos->s_star_down[1*N*M+j*M+M-1] = infos->s_star_down[1*N*M+j*M+M-2];
        infos->s_star_down[2*N*M+j*M+M-1] = infos->s_star_down[2*N*M+j*M+M-2];
        infos->s_star_down[3*N*M+j*M+M-1] = infos->s_star_down[3*N*M+j*M+M-2];
    }

    //  walls
    for (int i = 1; i < M-1; i++){
        infos->s_star_up[0*N*M+(N-1)*M+i] = infos->s_star_up[0*N*M+(N-2)*M+i];
        infos->s_star_up[1*N*M+(N-1)*M+i] = infos->s_star_up[1*N*M+(N-2)*M+i];
        infos->s_star_up[2*N*M+(N-1)*M+i] = 0.;
        infos->s_star_up[3*N*M+(N-1)*M+i] = infos->s_star_up[3*N*M+(N-2)*M+i];

        infos->s_star_down[0*N*M+0*M+i] = infos->s_star_down[0*N*M+1*M+i];
        infos->s_star_down[1*N*M+0*M+i] = infos->s_star_down[1*N*M+1*M+i];
        infos->s_star_down[2*N*M+0*M+i] = 0.;
        infos->s_star_down[3*N*M+0*M+i] = infos->s_star_down[3*N*M+1*M+i];
    }

    // frontière commune upstream
    for (int i = 0; i < infos->i_lead; i++){
        for (int k = 0; k < 4; k++)
            infos->s_star_up[k*N*M+0*M+i] = infos->s_star_down[k*N*M+(N-1)*M+i];
    }

    //airfoil
    for(int i = infos->i_lead; i < infos->i_trail; i++){
        // up : j=0
        double rho = infos->s_star_up[0*N*M+1*M+i];
        double u   = infos->s_star_up[1*N*M+1*M+i] / rho;
        double v   = infos->s_star_up[2*N*M+1*M+i] / rho;
        
        // normale à la paroi (perpendiculaire à la tangente (dx, dy_updi))
        double dx = infos->dxdi[i];
        double dy = infos->dydi_up[0*M+i];
        double norm2 = dx*dx + dy*dy;

        double c_y0 = (u*dx + v*dy)/norm2;

        infos->s_star_up[0*N*M+0*M+i] = rho;
        infos->s_star_up[1*N*M+0*M+i] = rho *c_y0*dx ;
        infos->s_star_up[2*N*M+0*M+i] = rho *c_y0*dy ;
        infos->s_star_up[3*N*M+0*M+i] = infos->s_star_up[3*N*M+1*M+i];

        // down : j=N-1, voisin = j=N-2
        rho = infos->s_star_down[0*N*M+(N-2)*M+i];
        u   = infos->s_star_down[1*N*M+(N-2)*M+i] / rho;
        v   = infos->s_star_down[2*N*M+(N-2)*M+i] / rho;
        
        dy = infos->dydi_down[(N-1)*M+i];
        norm2 = dx*dx + dy*dy;
        c_y0 = (u*dx + v*dy)/norm2;
        
        infos->s_star_down[0*N*M+(N-1)*M+i] = rho;
        infos->s_star_down[1*N*M+(N-1)*M+i] = rho *c_y0*dx ;
        infos->s_star_down[2*N*M+(N-1)*M+i] = rho *c_y0*dy ;
        infos->s_star_down[3*N*M+(N-1)*M+i] = infos->s_star_down[3*N*M+(N-2)*M+i];
    }

    // frontière commune downstream
    for (int i = infos->i_trail+1; i < M; i++){
        for (int k = 0; k < 4; k++)
            infos->s_star_up[k*N*M+0*M+i] = infos->s_star_down[k*N*M+(N-1)*M+i];
    }
}


// maccormac resolution

double terme_diffusion_fct(struct grid_and_data *infos, double*s_a, double*s_b,int k,int j,int i, int j_plus_1, int j_minus_1){

    // I use s_a == s_b for every cases except the ones at the common border, where s_a = s_down and s_b = s_up = related to j+1
    int N = infos->N;
    int M = infos->M;
    double eps_min = infos->eps_min;
    double C = infos->C;
    double c = infos->c;

    int i_plus_1 = i+1;
    int i_minus_1 = i-1;

    double delta_u_i_plus_1  = s_a[1*M*N+j*M+i_plus_1]/ s_a[0*M*N+j*M+i_plus_1] - s_a[1*M*N+j*M+i]/ s_a[0*M*N+j*M+i];
    double delta_u_i_minus_1 = s_a[1*M*N+j*M+i]/ s_a[0*M*N+j*M+i] - s_a[1*M*N+j*M+i_minus_1]/ s_a[0*M*N+j*M+i_minus_1];
    double delta_u_j_plus_1  = s_b[1*M*N+j_plus_1*M+i]/ s_b[0*M*N+j_plus_1*M+i] - s_a[1*M*N+j*M+i]/ s_a[0*M*N+j*M+i];
    double delta_u_j_minus_1 = s_a[1*M*N+j*M+i]/ s_a[0*M*N+j*M+i] - s_a[1*M*N+j_minus_1*M+i]/ s_a[0*M*N+j_minus_1*M+i];

    double delta_v_i_plus_1  = s_a[2*M*N+j*M+i_plus_1]/ s_a[0*M*N+j*M+i_plus_1] - s_a[2*M*N+j*M+i]/ s_a[0*M*N+j*M+i];
    double delta_v_i_minus_1 = s_a[2*M*N+j*M+i]/ s_a[0*M*N+j*M+i] - s_a[2*M*N+j*M+i_minus_1]/ s_a[0*M*N+j*M+i_minus_1];
    double delta_v_j_plus_1  = s_b[2*M*N+j_plus_1*M+i]/ s_b[0*M*N+j_plus_1*M+i] - s_a[2*M*N+j*M+i]/ s_a[0*M*N+j*M+i];
    double delta_v_j_minus_1 = s_a[2*M*N+j*M+i]/ s_a[0*M*N+j*M+i] - s_a[2*M*N+j_minus_1*M+i]/ s_a[0*M*N+j_minus_1*M+i];

    double eps_i_plus_half  = eps_min + C*c*sqrt(delta_u_i_plus_1*delta_u_i_plus_1 + delta_v_i_plus_1*delta_v_i_plus_1);
    double eps_i_minus_half = eps_min + C*c*sqrt(delta_u_i_minus_1*delta_u_i_minus_1 + delta_v_i_minus_1*delta_v_i_minus_1);

    double eps_j_plus_half  = eps_min + C*c*sqrt(delta_u_j_plus_1*delta_u_j_plus_1 + delta_v_j_plus_1*delta_v_j_plus_1);
    double eps_j_minus_half = eps_min + C*c*sqrt(delta_u_j_minus_1*delta_u_j_minus_1 + delta_v_j_minus_1*delta_v_j_minus_1);

    double terme_diffusion_i = eps_i_plus_half*(s_a[k*M*N+j*M+i_plus_1]-s_a[k*M*N+j*M+i]) -  eps_i_minus_half*(s_a[k*M*N+j*M+i]-s_a[k*M*N+j*M+i_minus_1]);
    double terme_diffusion_j = eps_j_plus_half*(s_b[k*M*N+j_plus_1*M+i]-s_a[k*M*N+j*M+i]) -  eps_j_minus_half*(s_a[k*M*N+j*M+i]-s_a[k*M*N+j_minus_1*M+i]);

    double terme_diffusion = terme_diffusion_i + terme_diffusion_j;

    return terme_diffusion;

}

void mac_cormack_forward(struct grid_and_data *infos){
    int N = infos->N;
    int M = infos->M;
    double dFidi;
    double dFjdj;
    
    double eps_min = infos->eps_min;
    double terme_diffusion;
    double terme_tot;
    //init FS

    for (int j = 0; j < N; j++){
        for (int i = 0; i < M; i++){

            double rho_up = infos->s_up[0*M*N+j*M+i];
            double u_up   = infos->s_up[1*M*N+j*M+i]/rho_up;
            double v_up   = infos->s_up[2*M*N+j*M+i]/rho_up;
            double U0_up  = infos->s_up[3*M*N+j*M+i]/rho_up;
            double p_up    = (infos->gamma - 1.) * rho_up * (U0_up - 0.5*(u_up*u_up + v_up*v_up));
            double H0_up = U0_up + p_up/rho_up;

            double rho_down = infos->s_down[0*M*N+j*M+i];
            double u_down   = infos->s_down[1*M*N+j*M+i]/rho_down;
            double v_down   = infos->s_down[2*M*N+j*M+i]/rho_down;
            double U0_down  = infos->s_down[3*M*N+j*M+i]/rho_down;
            double p_down   = (infos->gamma - 1.) * rho_down * (U0_down - 0.5*(u_down*u_down + v_down*v_down));
            double H0_down = U0_down + p_down/rho_down;

            double F_x_up_0 = infos->s_up[1*M*N+j*M+i];
            double F_x_up_1 = infos->s_up[1*M*N+j*M+i]*u_up + p_up;
            double F_x_up_2 = infos->s_up[1*M*N+j*M+i]*v_up;
            double F_x_up_3  = infos->s_up[1*M*N+j*M+i]*H0_up;

            double F_x_down_0 = infos->s_down[1*M*N+j*M+i];
            double F_x_down_1 = infos->s_down[1*M*N+j*M+i]*u_down + p_down;
            double F_x_down_2 = infos->s_down[1*M*N+j*M+i]*v_down;
            double F_x_down_3 = infos->s_down[1*M*N+j*M+i]*H0_down;

            double F_y_up_0 = infos->s_up[2*M*N+j*M+i];
            double F_y_up_1 = infos->s_up[2*M*N+j*M+i]*u_up;
            double F_y_up_2 = infos->s_up[2*M*N+j*M+i]*v_up + p_up;
            double F_y_up_3 = infos->s_up[2*M*N+j*M+i]*H0_up;

            double F_y_down_0 = infos->s_down[2*M*N+j*M+i];
            double F_y_down_1 = infos->s_down[2*M*N+j*M+i]*u_down;
            double F_y_down_2 = infos->s_down[2*M*N+j*M+i]*v_down + p_down;
            double F_y_down_3 = infos->s_down[2*M*N+j*M+i]*H0_down;

            
            double F_x_up_arr[4]   = {F_x_up_0,   F_x_up_1,   F_x_up_2,   F_x_up_3};
            double F_x_down_arr[4] = {F_x_down_0, F_x_down_1, F_x_down_2, F_x_down_3};
            double F_y_up_arr[4]   = {F_y_up_0,   F_y_up_1,   F_y_up_2,   F_y_up_3};
            double F_y_down_arr[4] = {F_y_down_0, F_y_down_1, F_y_down_2, F_y_down_3};

            for (int k = 0; k < 4; k++){
                infos->F_i_up[k*M*N + j*M + i]   = infos->dydj_up[j*M + i]   * F_x_up_arr[k];
                infos->F_i_down[k*M*N + j*M + i] = infos->dydj_down[j*M + i] * F_x_down_arr[k];

                infos->F_j_up[k*M*N + j*M + i]   =
                    -infos->dydt_up[j*M + i]*infos->dxdi[i]*infos->s_up[k*M*N + j*M + i]
                    - infos->dydi_up[j*M + i]*F_x_up_arr[k]
                    + infos->dxdi[i]*F_y_up_arr[k];

                infos->F_j_down[k*M*N + j*M + i] =
                    -infos->dydt_down[j*M + i]*infos->dxdi[i]*infos->s_down[k*M*N + j*M + i]
                    - infos->dydi_down[j*M + i]*F_x_down_arr[k]
                    + infos->dxdi[i]*F_y_down_arr[k];
            }
        }
    }

    for (int i = 0; i<infos->i_lead;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_up[k*N*M+0*M+i] = infos->F_i_down[k*N*M+(N-1)*M+i];
            infos->F_j_up[k*N*M+0*M+i] = infos->F_j_down[k*N*M+(N-1)*M+i];            
        }
    }

    for (int i = infos->i_trail+1; i<M;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_up[k*N*M+0*M+i] = infos->F_i_down[k*N*M+(N-1)*M+i];
            infos->F_j_up[k*N*M+0*M+i] = infos->F_j_down[k*N*M+(N-1)*M+i];     
        }
    }

    //Mac cormac




    //down part wall, 


    for (int i = 1; i <M-1;i++ ){
        int j = 0;
        for (int k = 0; k<4;k++){
            dFidi = infos->F_i_down[k*N*M+j*M+i+1]-infos->F_i_down[k*M*N+j*M+i]; 
            dFjdj = infos->F_j_down[k*M*N+j*M+i] - infos->F_j_down[k*N*M+(j+1)*M+i];  // ← inversé

            // j-1 devient j+1
            //terme_diffusion = eps_min * (infos->s_down[k*N*M+j*M+i+1] - 2*infos->s_down[k*M*N+j*M+i] + infos->s_down[k*N*M+j*M+i-1]); // j-1 devient j+1
            terme_diffusion = terme_diffusion_fct(infos, infos->s_down,infos->s_down,k, j, i, j+1, j+1);
             
            terme_tot = -dFidi + terme_diffusion;
            infos->s_star_down[k*M*N+j*M+i] = infos->s_down[k*M*N+j*M+i]*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + infos->dt*terme_tot/infos->D_down_t_plus_1[j*M+i];        
        }
    }
    // down part
    for (int j = 1; j < N-1; j++){
        for (int i = 1; i < M-1;i++){
            for (int k = 0; k<4;k++){
                dFidi = infos->F_i_down[k*N*M+j*M+i+1]-infos->F_i_down[k*M*N+j*M+i];
                dFjdj = infos->F_j_down[k*N*M+(j+1)*M+i]-infos->F_j_down[k*M*N+j*M+i];

                //terme_diffusion = eps_min * (infos->s_down[k*N*M+j*M+i+1] - 2*infos->s_down[k*M*N+j*M+i] + infos->s_down[k*N*M+j*M+i-1])+ eps_min * (infos->s_down[k*N*M+(j+1)*M+i] - 2*infos->s_down[k*M*N+j*M+i] + infos->s_down[k*N*M+(j-1)*M+i]); 
                terme_diffusion = terme_diffusion_fct(infos, infos->s_down,infos->s_down,k, j, i, j+1, j-1);
                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_star_down[k*M*N+j*M+i] = infos->s_down[k*M*N+j*M+i]*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + infos->dt*terme_tot/infos->D_down_t_plus_1[j*M+i];
                
            }
        }
    }


    // frontière upstream  
    for (int i = 1; i < infos->i_lead;i++){
        for (int k = 0; k < 4; k++){
            int j = N-1;

            dFidi = infos->F_i_down[k*N*M+j*M+i+1]-infos->F_i_down[k*M*N+j*M+i];
            dFjdj = infos->F_j_up[k*N*M+1*M+i]-infos->F_j_down[k*M*N+j*M+i];
            //dFjdj = 0.0;  // TEST
            
            terme_diffusion = terme_diffusion_fct(infos,infos->s_down, infos->s_up, k, j, i, 1, j-1 );
            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            infos->s_star_down[k*M*N+j*M+i] = infos->s_down[k*M*N+j*M+i]*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + infos->dt*terme_tot/infos->D_down_t_plus_1[j*M+i];

            infos->F_j_up[k*N*M+0*M+i] = infos->F_j_down[k*M*N+j*M+i];
            infos->F_i_up[k*N*M+0*M+i] = infos->F_i_down[k*M*N+j*M+i];
            infos->s_star_up[k*N*M+0*M+i] = infos->s_star_down[k*M*N+j*M+i];


        }
    }


    //up and downwing:

    for (int i = infos->i_lead; i <= infos->i_trail;i++){
        for (int k = 0; k < 4; k++){

            int j = N-1;
            //down wing, j+1 devient j-1

            dFidi = infos->F_i_down[k*N*M+j*M+i+1]-infos->F_i_down[k*M*N+j*M+i];


            terme_diffusion = terme_diffusion_fct(infos, infos->s_down,infos->s_down, k, j, i, j-1, j-1);
             
            terme_tot = -dFidi  + terme_diffusion;
            infos->s_star_down[k*N*M+(N-1)*M+i] = infos->s_down[k*N*M+(N-1)*M+i]*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + infos->dt*terme_tot/infos->D_down_t_plus_1[(N-1)*M+i];

            //upwing j-1 devient j+1
            j = 0;

            dFidi = infos->F_i_up[k*N*M+j*M+i+1]-infos->F_i_up[k*M*N+j*M+i];

            terme_diffusion = terme_diffusion_fct(infos, infos->s_up,infos->s_up, k, j, i, j+1, j+1);
            terme_tot = -dFidi + terme_diffusion;
            //printf("forw : i = %i, terme tot down = %f", i , terme_tot);

            infos->s_star_up[k*M*N+j*M+i] = infos->s_up[k*M*N+j*M+i]*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i] + infos->dt*terme_tot/infos->D_up_t_plus_1[j*M+i];

        }
    }

    //frontière downstream; 
    for (int i = infos->i_trail+1; i < M-1;i++){
        for (int k = 0; k < 4; k++){
            int j = N-1;
            dFidi = infos->F_i_down[k*N*M+(N-1)*M+i+1]-infos->F_i_down[k*N*M+(N-1)*M+i];
            dFjdj = infos->F_j_up[k*N*M+1*M+i]-infos->F_j_down[k*N*M+(N-1)*M+i];

            terme_diffusion = terme_diffusion_fct(infos,infos->s_down, infos->s_up, k, j, i, 1, j-1 );
            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            
            infos->s_star_down[k*N*M+(N-1)*M+i] = infos->s_down[k*N*M+(N-1)*M+i]*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + infos->dt*terme_tot/infos->D_down_t_plus_1[(N-1)*M+i];

            infos->F_j_up[k*N*M+0*M+i] = infos->F_j_down[k*N*M+(N-1)*M+i];
            infos->F_i_up[k*N*M+0*M+i] = infos->F_i_down[k*N*M+(N-1)*M+i];
            infos->s_star_up[k*N*M+0*M+i] = infos->s_star_down[k*N*M+(N-1)*M+i];
        }

    }


    //upper part exept the wall

    for (int j = 1; j < N-1; j++){
        for (int i = 1; i < M-1;i++){
            for (int k = 0; k<4;k++){

                dFidi = infos->F_i_up[k*N*M+j*M+i+1]-infos->F_i_up[k*M*N+j*M+i];
                dFjdj = infos->F_j_up[k*N*M+(j+1)*M+i]-infos->F_j_up[k*M*N+j*M+i];

                //terme_diffusion = eps_min * (infos->s_up[k*N*M+j*M+i+1] - 2*infos->s_up[k*M*N+j*M+i] + infos->s_up[k*N*M+j*M+i-1])+ eps_min * (infos->s_up[k*N*M+(j+1)*M+i] - 2*infos->s_up[k*M*N+j*M+i] + infos->s_up[k*N*M+(j-1)*M+i]);
                terme_diffusion = terme_diffusion_fct(infos, infos->s_up,infos->s_up, k, j, i, j+1, j-1);

                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_star_up[k*M*N+j*M+i] = infos->s_up[k*M*N+j*M+i]*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i] + infos->dt*terme_tot/infos->D_up_t_plus_1[j*M+i];

            }
        }
    }

    //upper wall, j+1 devient j-1
    
    for (int i = 1; i < M-1;i++){
        int j = N-1;
        for (int k = 0; k<4;k++){
            dFidi = infos->F_i_up[k*N*M+j*M+i+1]-infos->F_i_up[k*M*N+j*M+i];

            //terme_diffusion = eps_min * (infos->s_up[k*N*M+j*M+i+1] - 2*infos->s_up[k*M*N+j*M+i] + infos->s_up[k*N*M+j*M+i-1]);
            terme_diffusion = terme_diffusion_fct(infos, infos->s_up,infos->s_up, k, j, i, j-1, j-1);

            terme_tot = -dFidi + terme_diffusion;
            infos->s_star_up[k*M*N+j*M+i] = infos->s_up[k*M*N+j*M+i]*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i] + infos->dt*terme_tot/infos->D_up_t_plus_1[j*M+i];
        }
    }


}

void mac_cormack_backward_save(struct grid_and_data *infos){
    int N = infos->N;
    int M = infos->M;
    double dFidi;
    double dFjdj;
    double terme_diffusion;
    double terme_tot;
    double epsilon;

    //init les FS star

    for (int j = 0; j < N; j++){
        for (int i = 0; i < M; i++){
            double rho_star_up = infos->s_star_up[0*M*N+j*M+i];
            double u_star_up   = infos->s_star_up[1*M*N+j*M+i]/rho_star_up;
            double v_star_up   = infos->s_star_up[2*M*N+j*M+i]/rho_star_up;
            double U0_star_up  = infos->s_star_up[3*M*N+j*M+i]/rho_star_up;
            double p_star_up   = (infos->gamma - 1.) * rho_star_up * (U0_star_up - 0.5*(u_star_up*u_star_up + v_star_up*v_star_up));
            double H0_star_up = U0_star_up + p_star_up/rho_star_up;

            double rho_star_down = infos->s_star_down[0*M*N+j*M+i];
            double u_star_down   = infos->s_star_down[1*M*N+j*M+i]/rho_star_down;
            double v_star_down   = infos->s_star_down[2*M*N+j*M+i]/rho_star_down;
            double U0_star_down  = infos->s_star_down[3*M*N+j*M+i]/rho_star_down;
            double p_star_down   = (infos->gamma - 1.) * rho_star_down * (U0_star_down - 0.5*(u_star_down*u_star_down + v_star_down*v_star_down));
            double H0_star_down = U0_star_down + p_star_down/rho_star_down;

            double F_x_star_up_0 = infos->s_star_up[1*M*N+j*M+i];
            double F_x_star_up_1 = infos->s_star_up[1*M*N+j*M+i]*u_star_up+  p_star_up;
            double F_x_star_up_2 = infos->s_star_up[1*M*N+j*M+i]*v_star_up;
            double F_x_star_up_3 = infos->s_star_up[1*M*N+j*M+i]*H0_star_up;

            double F_x_star_down_0 = infos->s_star_down[1*M*N+j*M+i];
            double F_x_star_down_1 = infos->s_star_down[1*M*N+j*M+i]*u_star_down +p_star_down;
            double F_x_star_down_2 = infos->s_star_down[1*M*N+j*M+i]*v_star_down;
            double F_x_star_down_3 = infos->s_star_down[1*M*N+j*M+i]*H0_star_down;

            double F_y_star_up_0 = infos->s_star_up[2*M*N+j*M+i];
            double F_y_star_up_1 = infos->s_star_up[2*M*N+j*M+i]*u_star_up;
            double F_y_star_up_2 = infos->s_star_up[2*M*N+j*M+i]*v_star_up +p_star_up;
            double F_y_star_up_3 = infos->s_star_up[2*M*N+j*M+i]*H0_star_up;

            double F_y_star_down_0 = infos->s_star_down[2*M*N+j*M+i];
            double F_y_star_down_1 = infos->s_star_down[2*M*N+j*M+i]*u_star_down;
            double F_y_star_down_2 = infos->s_star_down[2*M*N+j*M+i]*v_star_down+ p_star_down;
            double F_y_star_down_3 = infos->s_star_down[2*M*N+j*M+i]*H0_star_down;

                        
            double F_x_star_up_arr[4]   = {F_x_star_up_0,   F_x_star_up_1,   F_x_star_up_2,   F_x_star_up_3};
            double F_x_star_down_arr[4] = {F_x_star_down_0, F_x_star_down_1, F_x_star_down_2, F_x_star_down_3};
            double F_y_star_up_arr[4]   = {F_y_star_up_0,   F_y_star_up_1,   F_y_star_up_2,   F_y_star_up_3};
            double F_y_star_down_arr[4] = {F_y_star_down_0, F_y_star_down_1, F_y_star_down_2, F_y_star_down_3};

            for (int k = 0; k < 4; k++){
                infos->F_i_star_up[k*M*N + j*M + i]   = infos->dydj_up[j*M + i]   * F_x_star_up_arr[k];
                infos->F_i_star_down[k*M*N + j*M + i] = infos->dydj_down[j*M + i] * F_x_star_down_arr[k];

                infos->F_j_star_up[k*M*N + j*M + i]   =
                    -infos->dydt_up[j*M + i]*infos->dxdi[i]*infos->s_star_up[k*M*N + j*M + i]
                    - infos->dydi_up[j*M + i]*F_x_star_up_arr[k]
                    + infos->dxdi[i]*F_y_star_up_arr[k];

                infos->F_j_star_down[k*M*N + j*M + i] =
                    -infos->dydt_down[j*M + i]*infos->dxdi[i]*infos->s_star_down[k*M*N + j*M + i]
                    - infos->dydi_down[j*M + i]*F_x_star_down_arr[k]
                    + infos->dxdi[i]*F_y_star_down_arr[k];
            }

        }
    }

    for (int i = 0; i<infos->i_lead;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_star_up[k*N*M+0*M+i] = infos->F_i_star_down[k*N*M+(N-1)*M+i];
            infos->F_j_star_up[k*N*M+0*M+i] = infos->F_j_star_down[k*N*M+(N-1)*M+i];            
        }
    }

    for (int i = infos->i_trail+1; i<M;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_star_up[k*N*M+0*M+i] = infos->F_i_star_down[k*N*M+(N-1)*M+i];
            infos->F_j_star_up[k*N*M+0*M+i] = infos->F_j_star_down[k*N*M+(N-1)*M+i];    
        }
    }


    // down wall j-1 devient j+1

    for (int i = 1; i < M-1;i++){
        for (int k = 0; k<4;k++){
            int j = 0;
            dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];
            //dFjdj = infos->F_j_star_down[k*N*M+(j+1)*M+i] - infos->F_j_star_down[k*M*N+j*M+i];  


            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_down[k*N*M+j*M+i+1] - 2*infos->s_star_down[k*M*N+j*M+i] + infos->s_star_down[k*N*M+j*M+i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j+1, j+1);

            terme_tot = -dFidi + terme_diffusion;
            infos->s_down[k*N*M+0*M+i] = 0.5*(infos->s_star_down[k*N*M+0*M+i] + infos->s_down[k*N*M+0*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i]+ 0.5*infos->dt/infos->D_down_t_plus_1[0*M+i] * terme_tot;
        }
    }


    //down part
    for (int j = 1; j < N-1; j++){
        for (int i = 1; i < M-1;i++){
            for (int k = 0; k<4;k++){

                dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];
                dFjdj = infos->F_j_star_down[k*M*N+j*M+i]-infos->F_j_star_down[k*N*M+(j-1)*M+i];
                //dFjdj = infos->F_j_star_down[k*N*M+(j+1)*M+i] - infos->F_j_star_down[k*M*N+j*M+i];  // ← inversé

                epsilon = infos->eps_min;
                //terme_diffusion = epsilon * (infos->s_star_down[k*N*M+j*M+i+1] - 2*infos->s_star_down[k*M*N+j*M+i] + infos->s_star_down[k*N*M+j*M+i-1])+ epsilon * (infos->s_star_down[k*N*M+(j-1)*M+i] - 2*infos->s_star_down[k*M*N+j*M+i] + infos->s_star_down[k*N*M+(j+1)*M+i]);
                
                terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j+1, j-1);
                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_down[k*M*N+j*M+i] = 0.5*(infos->s_star_down[k*M*N+j*M+i] + infos->s_down[k*M*N+j*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i]+ 0.5*infos->dt/infos->D_down_t_plus_1[j*M+i] * terme_tot;

            }
        }
    }

    // frontière upstream 
    for (int i = 1; i < infos->i_lead;i++){
        int j = N-1;
        for (int k = 0; k < 4; k++){
            dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];
            dFjdj = infos->F_j_star_down[k*M*N+j*M+i]-infos->F_j_star_up[k*N*M+1*M+i];
            //dFjdj = 0.0;  // TEST

            epsilon = infos->eps_min;
            terme_diffusion  = terme_diffusion_fct(infos,infos->s_star_down, infos->s_star_up, k, j, i, 1, j-1 );

            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            infos->s_down[k*M*N+j*M+i] = 0.5*(infos->s_star_down[k*M*N+j*M+i] + infos->s_down[k*M*N+j*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_down_t_plus_1[j*M+i] * terme_tot;

            infos->F_j_star_up[k*N*M+0*M+i] = infos->F_j_star_down[k*M*N+j*M+i];
            infos->F_i_star_up[k*N*M+0*M+i] = infos->F_i_star_down[k*M*N+j*M+i];
            //infos->s_up[k*N*M+0*M+i] = infos->s_down[k*N*M+(N-1)*M+i];
        }
    }

    //up and down wing:
    
    for (int i = infos->i_lead; i < infos->i_trail+1;i++){
        for (int k = 0; k < 4; k++){
            int j = N-1;
            //down wing, j+1 devient j-1

            dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];

            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_down[k*N*M+j*M+i+1] - 2*infos->s_star_down[k*M*N+j*M+i] + infos->s_star_down[k*N*M+j*M+i-1]);
            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j-1, j-1);

            terme_tot = -dFidi + terme_diffusion;
            infos->s_down[k*M*N+j*M+i] = 0.5*(infos->s_star_down[k*M*N+j*M+i] + infos->s_down[k*M*N+j*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_down_t_plus_1[j*M+i] * terme_tot;

            //upwing j-1 devient j+1
            j = 0;
            dFidi = infos->F_i_star_up[k*M*N+j*M+i]-infos->F_i_star_up[k*N*M+j*M+i-1];

            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_up[k*N*M+j*M+i+1] - 2*infos->s_star_up[k*M*N+j*M+i] + infos->s_star_up[k*N*M+j*M+i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,  infos->s_star_up, k, j, i, j+1, j+1);
            terme_tot = -dFidi + terme_diffusion;

            infos->s_up[k*M*N+j*M+i] = 0.5*(infos->s_star_up[k*M*N+j*M+i] + infos->s_up[k*M*N+j*M+i])*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i]+ 0.5*infos->dt/infos->D_up_t_plus_1[j*M+i] * terme_tot;

        }

    }
    //frontière downstream;
    for (int i = infos->i_trail+1; i < M-1;i++){
        for (int k = 0; k < 4; k++){
            int j = N-1;
            dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];
            dFjdj = infos->F_j_star_down[k*M*N+j*M+i]-infos->F_j_star_up[k*N*M+1*M+i];
            //dFjdj = 0.0;  // TEST

            epsilon = infos->eps_min;
            terme_diffusion = terme_diffusion_fct(infos,infos->s_star_down, infos->s_star_up, k, j, i, 1, j-1 );

            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            infos->s_down[k*M*N+j*M+i] = 0.5*(infos->s_star_down[k*M*N+j*M+i] + infos->s_down[k*M*N+j*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_down_t_plus_1[j*M+i] * terme_tot;

            infos->F_j_star_up[k*N*M+0*M+i] = infos->F_j_star_down[k*M*N+j*M+i];
            infos->F_i_star_up[k*N*M+0*M+i] = infos->F_i_star_down[k*M*N+j*M+i];
            //infos->s_up[k*N*M+0*M+i] = infos->s_down[k*M*N+j*M+i];

        }
    }
    //upper part
    for (int j = 1; j < N-1; j++){
        for (int i = 1; i < M-1;i++){
            for (int k = 0; k<4;k++){

                dFidi = infos->F_i_star_up[k*M*N+j*M+i]-infos->F_i_star_up[k*N*M+j*M+i-1];
                dFjdj = infos->F_j_star_up[k*M*N+j*M+i]-infos->F_j_star_up[k*N*M+(j-1)*M+i];

                epsilon = infos->eps_min;
                //terme_diffusion = epsilon * (infos->s_star_up[k*N*M+j*M+i+1] - 2*infos->s_star_up[k*M*N+j*M+i] + infos->s_star_up[k*N*M+j*M+i-1])+ epsilon * (infos->s_star_up[k*N*M+(j+1)*M+i] - 2*infos->s_star_up[k*M*N+j*M+i] + infos->s_star_up[k*N*M+(j-1)*M+i]);
                terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,infos->s_star_up, k, j, i, j+1, j-1);
                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_up[k*M*N+j*M+i] = 0.5*(infos->s_star_up[k*M*N+j*M+i] + infos->s_up[k*M*N+j*M+i])*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_up_t_plus_1[j*M+i] * terme_tot;
            }
        }

    }

    // upper wall j+1 devient j-1

    for (int i = 1; i < M-1;i++){
        int j = N-1;
        for (int k = 0; k<4;k++){

            dFidi = infos->F_i_star_up[k*M*N+j*M+i]-infos->F_i_star_up[k*N*M+j*M+i-1];

            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_up[k*N*M+j*M+i+1] - 2*infos->s_star_up[k*M*N+j*M+i] + infos->s_star_up[k*N*M+j*M+i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,infos->s_star_up, k, j, i, j-1, j-1);
        
            terme_tot = -dFidi  + terme_diffusion;
            infos->s_up[k*M*N+j*M+i] = 0.5*(infos->s_star_up[k*M*N+j*M+i] + infos->s_up[k*M*N+j*M+i])*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_up_t_plus_1[j*M+i] * terme_tot;
        }
    }
}

void mac_cormack_backward(struct grid_and_data *infos){
    int N = infos->N;
    int M = infos->M;
    double dFidi;
    double dFjdj;
    double terme_diffusion;
    double terme_tot;
    double epsilon;

    //init les FS star

    for (int j = 0; j < N; j++){
        for (int i = 0; i < M; i++){
            double rho_star_up = infos->s_star_up[0*M*N+j*M+i];
            double u_star_up   = infos->s_star_up[1*M*N+j*M+i]/rho_star_up;
            double v_star_up   = infos->s_star_up[2*M*N+j*M+i]/rho_star_up;
            double U0_star_up  = infos->s_star_up[3*M*N+j*M+i]/rho_star_up;
            double p_star_up   = (infos->gamma - 1.) * rho_star_up * (U0_star_up - 0.5*(u_star_up*u_star_up + v_star_up*v_star_up));
            double H0_star_up = U0_star_up + p_star_up/rho_star_up;

            double rho_star_down = infos->s_star_down[0*M*N+j*M+i];
            double u_star_down   = infos->s_star_down[1*M*N+j*M+i]/rho_star_down;
            double v_star_down   = infos->s_star_down[2*M*N+j*M+i]/rho_star_down;
            double U0_star_down  = infos->s_star_down[3*M*N+j*M+i]/rho_star_down;
            double p_star_down   = (infos->gamma - 1.) * rho_star_down * (U0_star_down - 0.5*(u_star_down*u_star_down + v_star_down*v_star_down));
            double H0_star_down = U0_star_down + p_star_down/rho_star_down;

            double F_x_star_up_0 = infos->s_star_up[1*M*N+j*M+i];
            double F_x_star_up_1 = infos->s_star_up[1*M*N+j*M+i]*u_star_up+  p_star_up;
            double F_x_star_up_2 = infos->s_star_up[1*M*N+j*M+i]*v_star_up;
            double F_x_star_up_3 = infos->s_star_up[1*M*N+j*M+i]*H0_star_up;

            double F_x_star_down_0 = infos->s_star_down[1*M*N+j*M+i];
            double F_x_star_down_1 = infos->s_star_down[1*M*N+j*M+i]*u_star_down +p_star_down;
            double F_x_star_down_2 = infos->s_star_down[1*M*N+j*M+i]*v_star_down;
            double F_x_star_down_3 = infos->s_star_down[1*M*N+j*M+i]*H0_star_down;

            double F_y_star_up_0 = infos->s_star_up[2*M*N+j*M+i];
            double F_y_star_up_1 = infos->s_star_up[2*M*N+j*M+i]*u_star_up;
            double F_y_star_up_2 = infos->s_star_up[2*M*N+j*M+i]*v_star_up +p_star_up;
            double F_y_star_up_3 = infos->s_star_up[2*M*N+j*M+i]*H0_star_up;

            double F_y_star_down_0 = infos->s_star_down[2*M*N+j*M+i];
            double F_y_star_down_1 = infos->s_star_down[2*M*N+j*M+i]*u_star_down;
            double F_y_star_down_2 = infos->s_star_down[2*M*N+j*M+i]*v_star_down+ p_star_down;
            double F_y_star_down_3 = infos->s_star_down[2*M*N+j*M+i]*H0_star_down;

                        
            double F_x_star_up_arr[4]   = {F_x_star_up_0,   F_x_star_up_1,   F_x_star_up_2,   F_x_star_up_3};
            double F_x_star_down_arr[4] = {F_x_star_down_0, F_x_star_down_1, F_x_star_down_2, F_x_star_down_3};
            double F_y_star_up_arr[4]   = {F_y_star_up_0,   F_y_star_up_1,   F_y_star_up_2,   F_y_star_up_3};
            double F_y_star_down_arr[4] = {F_y_star_down_0, F_y_star_down_1, F_y_star_down_2, F_y_star_down_3};

            for (int k = 0; k < 4; k++){
                infos->F_i_star_up[k*M*N + j*M + i]   = infos->dydj_up[j*M + i]   * F_x_star_up_arr[k];
                infos->F_i_star_down[k*M*N + j*M + i] = infos->dydj_down[j*M + i] * F_x_star_down_arr[k];

                infos->F_j_star_up[k*M*N + j*M + i]   =
                    -infos->dydt_up[j*M + i]*infos->dxdi[i]*infos->s_star_up[k*M*N + j*M + i]
                    - infos->dydi_up[j*M + i]*F_x_star_up_arr[k]
                    + infos->dxdi[i]*F_y_star_up_arr[k];

                infos->F_j_star_down[k*M*N + j*M + i] =
                    -infos->dydt_down[j*M + i]*infos->dxdi[i]*infos->s_star_down[k*M*N + j*M + i]
                    - infos->dydi_down[j*M + i]*F_x_star_down_arr[k]
                    + infos->dxdi[i]*F_y_star_down_arr[k];
            }

        }
    }

    for (int i = 0; i<infos->i_lead;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_star_up[k*N*M+0*M+i] = infos->F_i_star_down[k*N*M+(N-1)*M+i];
            infos->F_j_star_up[k*N*M+0*M+i] = infos->F_j_star_down[k*N*M+(N-1)*M+i];            
        }
    }

    for (int i = infos->i_trail+1; i<M;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_star_up[k*N*M+0*M+i] = infos->F_i_star_down[k*N*M+(N-1)*M+i];
            infos->F_j_star_up[k*N*M+0*M+i] = infos->F_j_star_down[k*N*M+(N-1)*M+i];    
        }
    }


    // down wall j-1 devient j+1

    for (int i = 1; i < M-1;i++){
        for (int k = 0; k<4;k++){
            int j = 0;
            dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];
            //dFjdj = infos->F_j_star_down[k*N*M+(j+1)*M+i] - infos->F_j_star_down[k*M*N+j*M+i];  


            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_down[k*N*M+j*M+i+1] - 2*infos->s_star_down[k*M*N+j*M+i] + infos->s_star_down[k*N*M+j*M+i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j+1, j+1);

            terme_tot = -dFidi + terme_diffusion;
            infos->s_down[k*N*M+0*M+i] = 0.5*(infos->s_star_down[k*N*M+0*M+i] + infos->s_down[k*N*M+0*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i]+ 0.5*infos->dt/infos->D_down_t_plus_1[0*M+i] * terme_tot;
        }
    }


    //down part
    for (int j = 1; j < N-1; j++){
        for (int i = 1; i < M-1;i++){
            for (int k = 0; k<4;k++){

                dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];
                dFjdj = infos->F_j_star_down[k*M*N+j*M+i]-infos->F_j_star_down[k*N*M+(j-1)*M+i];
                //dFjdj = infos->F_j_star_down[k*N*M+(j+1)*M+i] - infos->F_j_star_down[k*M*N+j*M+i];  // ← inversé

                epsilon = infos->eps_min;
                //terme_diffusion = epsilon * (infos->s_star_down[k*N*M+j*M+i+1] - 2*infos->s_star_down[k*M*N+j*M+i] + infos->s_star_down[k*N*M+j*M+i-1])+ epsilon * (infos->s_star_down[k*N*M+(j-1)*M+i] - 2*infos->s_star_down[k*M*N+j*M+i] + infos->s_star_down[k*N*M+(j+1)*M+i]);
                
                terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j+1, j-1);
                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_down[k*M*N+j*M+i] = 0.5*(infos->s_star_down[k*M*N+j*M+i] + infos->s_down[k*M*N+j*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i]+ 0.5*infos->dt/infos->D_down_t_plus_1[j*M+i] * terme_tot;

            }
        }
    }

    // frontière upstream 
    for (int i = 1; i < infos->i_lead;i++){
        int j = N-1;
        for (int k = 0; k < 4; k++){
            dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];
            dFjdj = infos->F_j_star_down[k*M*N+j*M+i]-infos->F_j_star_down[k*N*M+(j-1)*M+i];
            //dFjdj = 0.0;  // TEST

            epsilon = infos->eps_min;
            terme_diffusion  = terme_diffusion_fct(infos,infos->s_star_down, infos->s_star_up, k, j, i, 1, j-1 );

            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            infos->s_down[k*M*N+j*M+i] = 0.5*(infos->s_star_down[k*M*N+j*M+i] + infos->s_down[k*M*N+j*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_down_t_plus_1[j*M+i] * terme_tot;

            infos->F_j_star_up[k*N*M+0*M+i] = infos->F_j_star_down[k*M*N+j*M+i];
            infos->F_i_star_up[k*N*M+0*M+i] = infos->F_i_star_down[k*M*N+j*M+i];
            //infos->s_up[k*N*M+0*M+i] = infos->s_down[k*N*M+(N-1)*M+i];
        }
    }

    //up and down wing:
    
    for (int i = infos->i_lead; i < infos->i_trail+1;i++){
        for (int k = 0; k < 4; k++){
            int j = N-1;
            //down wing, j+1 devient j-1

            dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];

            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_down[k*N*M+j*M+i+1] - 2*infos->s_star_down[k*M*N+j*M+i] + infos->s_star_down[k*N*M+j*M+i-1]);
            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j-1, j-1);

            terme_tot = -dFidi + terme_diffusion;
            infos->s_down[k*M*N+j*M+i] = 0.5*(infos->s_star_down[k*M*N+j*M+i] + infos->s_down[k*M*N+j*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_down_t_plus_1[j*M+i] * terme_tot;

            //upwing j-1 devient j+1
            j = 0;
            dFidi = infos->F_i_star_up[k*M*N+j*M+i]-infos->F_i_star_up[k*N*M+j*M+i-1];

            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_up[k*N*M+j*M+i+1] - 2*infos->s_star_up[k*M*N+j*M+i] + infos->s_star_up[k*N*M+j*M+i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,  infos->s_star_up, k, j, i, j+1, j+1);
            terme_tot = -dFidi + terme_diffusion;

            infos->s_up[k*M*N+j*M+i] = 0.5*(infos->s_star_up[k*M*N+j*M+i] + infos->s_up[k*M*N+j*M+i])*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i]+ 0.5*infos->dt/infos->D_up_t_plus_1[j*M+i] * terme_tot;

        }

    }
    //frontière downstream;
    for (int i = infos->i_trail+1; i < M-1;i++){
        for (int k = 0; k < 4; k++){
            int j = N-1;
            dFidi = infos->F_i_star_down[k*M*N+j*M+i]-infos->F_i_star_down[k*N*M+j*M+i-1];
            dFjdj = infos->F_j_star_down[k*M*N+j*M+i]-infos->F_j_star_down[k*N*M+(j-1)*M+i];
            //dFjdj = 0.0;  // TEST

            epsilon = infos->eps_min;
            terme_diffusion = terme_diffusion_fct(infos,infos->s_star_down, infos->s_star_up, k, j, i, 1, j-1 );

            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            infos->s_down[k*M*N+j*M+i] = 0.5*(infos->s_star_down[k*M*N+j*M+i] + infos->s_down[k*M*N+j*M+i])*infos->D_down[j*M+i]/infos->D_down_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_down_t_plus_1[j*M+i] * terme_tot;

            infos->F_j_star_up[k*N*M+0*M+i] = infos->F_j_star_down[k*M*N+j*M+i];
            infos->F_i_star_up[k*N*M+0*M+i] = infos->F_i_star_down[k*M*N+j*M+i];
            //infos->s_up[k*N*M+0*M+i] = infos->s_down[k*M*N+j*M+i];

        }
    }
    //upper part
    for (int j = 1; j < N-1; j++){
        for (int i = 1; i < M-1;i++){
            for (int k = 0; k<4;k++){

                dFidi = infos->F_i_star_up[k*M*N+j*M+i]-infos->F_i_star_up[k*N*M+j*M+i-1];
                dFjdj = infos->F_j_star_up[k*M*N+j*M+i]-infos->F_j_star_up[k*N*M+(j-1)*M+i];

                epsilon = infos->eps_min;
                //terme_diffusion = epsilon * (infos->s_star_up[k*N*M+j*M+i+1] - 2*infos->s_star_up[k*M*N+j*M+i] + infos->s_star_up[k*N*M+j*M+i-1])+ epsilon * (infos->s_star_up[k*N*M+(j+1)*M+i] - 2*infos->s_star_up[k*M*N+j*M+i] + infos->s_star_up[k*N*M+(j-1)*M+i]);
                terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,infos->s_star_up, k, j, i, j+1, j-1);
                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_up[k*M*N+j*M+i] = 0.5*(infos->s_star_up[k*M*N+j*M+i] + infos->s_up[k*M*N+j*M+i])*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_up_t_plus_1[j*M+i] * terme_tot;
            }
        }

    }

    // upper wall j+1 devient j-1

    for (int i = 1; i < M-1;i++){
        int j = N-1;
        for (int k = 0; k<4;k++){

            dFidi = infos->F_i_star_up[k*M*N+j*M+i]-infos->F_i_star_up[k*N*M+j*M+i-1];

            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_up[k*N*M+j*M+i+1] - 2*infos->s_star_up[k*M*N+j*M+i] + infos->s_star_up[k*N*M+j*M+i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,infos->s_star_up, k, j, i, j-1, j-1);
        
            terme_tot = -dFidi  + terme_diffusion;
            infos->s_up[k*M*N+j*M+i] = 0.5*(infos->s_star_up[k*M*N+j*M+i] + infos->s_up[k*M*N+j*M+i])*infos->D_up[j*M+i]/infos->D_up_t_plus_1[j*M+i] + 0.5*infos->dt/infos->D_up_t_plus_1[j*M+i] * terme_tot;
        }
    }
}


//recap_function

void solver_all_steps(double M_inf, double grow_rate_M_inf, double t_end ,double number_of_wrotten_steps, double resolution, double amplitude_oscillation, double frequence_oscillation, double angle_of_attack, bool hdd_storage){

    int result = mkdir("Animations", 0777);
    result = mkdir("binaryFiles", 0777);

    time_t begin = time(0); 
    struct grid_and_data infos = init_grid_and_data(M_inf, grow_rate_M_inf, t_end, number_of_wrotten_steps,resolution,  amplitude_oscillation, frequence_oscillation,angle_of_attack, hdd_storage);
    bounds_s_function set_bounds_s =  set_bounds_s_subsonic  ;
    bounds_s_function set_bounds_s_star = set_bounds_s_star_subsonic;

    if (infos.M_inf > 1. || grow_rate_M_inf !=0.){
        set_bounds_s =  set_bounds_s_supersonic  ;
        set_bounds_s_star = set_bounds_s_star_supersonic;
    }
    set_x_grid(&infos);
    leading_and_trailing_edge_finder(&infos);   
    set_y_grid(&infos);
    set_init_s(&infos);
    set_bounds_s(&infos);

    
    char * filename = generate_filename(&infos);
    FILE * file_binary = fopen(filename, "wb"); 
    write_first_two_lines(&infos, file_binary);

    double steps = 0.;
    double steps_already_done = 0.;

    while(infos.t < t_end){
        mac_cormack_forward(&infos);
        set_bounds_s_star(&infos);
        mac_cormack_backward(&infos);
        set_bounds_s(&infos);

        if (steps > infos.write_interval_steps){
            write_one_time_step(&infos, file_binary);
            steps = 0.;
            //printf("%g/100 \n", 100.*infos.t/infos.t_end);
        }

        if (infos.mu != 0 ||  infos.AOA_grow_rate !=0.){
            infos.AOA_current += infos.AOA_grow_rate * infos.dt;
            set_y_grid(&infos);
        }

        if (infos.AOA_current > infos.AOA){
            infos.AOA_grow_rate =0.;
            infos.AOA_current  = infos.AOA;
        }

        infos.t += infos.dt;
        steps +=1.;
        steps_already_done+=1.;
        infos.M_inf += infos.grow_rate_M_inf * infos.dt; /// ????? marche pas
    }

    fprintf(file_binary, "end \n");
    fclose(file_binary);
    free_data(&infos);

    time_t end = time(0); 
    printf("Done with %s after %ld seconds with %f timesteps\n", filename,  end - begin, steps_already_done);

}

#endif 