//
//  solver.h
//  MECA2660 projet 2026


#ifndef solver_h
#define solver_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  
#include "grid_and_data.h"
#include "debug_functions.h"


void set_init_s(struct grid_and_data *infos){

    double gamma  = infos->gamma;
    double rho_inf = infos->rho_inf;
    double p_inf   = infos->p_inf;
    double u_inf   = infos->M_inf * infos->c_sound_inf;

    double U0_inf = p_inf / ((gamma - 1.) * rho_inf) + 0.5 * u_inf * u_inf;      

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            infos->s_up[0][j][i] = rho_inf;
            infos->s_up[1][j][i] = rho_inf * u_inf;
            infos->s_up[2][j][i] = 0.;
            infos->s_up[3][j][i] = rho_inf * U0_inf;

            infos->s_down[0][j][i] = rho_inf;
            infos->s_down[1][j][i] = rho_inf * u_inf;
            infos->s_down[2][j][i] = 0.;
            infos->s_down[3][j][i] = rho_inf * U0_inf;

            for (int k = 0; k < 4; k++){
                infos->s_star_up[k][j][i]   = infos->s_up[k][j][i];
                infos->s_star_down[k][j][i] = infos->s_down[k][j][i];
            }
        }
    }




}

void set_bounds_s(struct grid_and_data *infos ){

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


    for (int j = 0; j< infos->N; j++){
        //inlet conditions

        rho_0 = infos->s_up[0][j][1]; 
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        arg = u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.) * (p_inf/rho_inf - p_0/rho_0);
        u_0 = sqrt(fabs(arg));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);

        infos->s_up[0][j][0] = rho_0;
        infos->s_up[1][j][0] = rho_0*u_0;
        infos->s_up[2][j][0] = rho_0*v_0;
        infos->s_up[3][j][0] = rho_0*U_0;

        rho_0 = infos->s_down[0][j][1];
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        arg = u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.) * (p_inf/rho_inf - p_0/rho_0);
        u_0 = sqrt(fabs(arg));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);

        infos->s_down[0][j][0]= rho_0;
        infos->s_down[1][j][0]= rho_0*u_0;
        infos->s_down[2][j][0]= rho_0*v_0;
        infos->s_down[3][j][0]= rho_0*U_0;

        rho_end = infos->s_up[0][j][infos->M-2];
        u_end = infos->s_up[1][j][infos->M-2]/rho_end;
        v_end = infos->s_up[2][j][infos->M-2]/rho_end;

        U0_end = p_inf/((infos->gamma-1.)*rho_end) + 0.5*(u_end*u_end+v_end*v_end);

        infos->s_up[0][j][infos->M-1] = infos->s_up[0][j][infos->M-2];
        infos->s_up[1][j][infos->M-1] = infos->s_up[1][j][infos->M-2];
        infos->s_up[2][j][infos->M-1] = infos->s_up[2][j][infos->M-2];
        infos->s_up[3][j][infos->M-1] = rho_end*U0_end;

        rho_end = infos->s_down[0][j][infos->M-2];
        u_end = infos->s_down[1][j][infos->M-2]/rho_end;
        v_end = infos->s_down[2][j][infos->M-2]/rho_end;
        U0_end = p_inf/((infos->gamma-1.)*rho_end) + 0.5*(u_end*u_end+v_end*v_end);

        infos->s_down[0][j][infos->M-1] = infos->s_down[0][j][infos->M-2];
        infos->s_down[1][j][infos->M-1] = infos->s_down[1][j][infos->M-2];
        infos->s_down[2][j][infos->M-1] = infos->s_down[2][j][infos->M-2];
        infos->s_down[3][j][infos->M-1] = rho_end*U0_end;
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

    for (int i = 1; i <infos->M-1; i++){

        infos->s_up[0][infos->N-1][i] = infos->s_up[0][infos->N-2][i];
        infos->s_up[1][infos->N-1][i] = infos->s_up[1][infos->N-2][i];
        infos->s_up[2][infos->N-1][i] = 0.0;
        infos->s_up[3][infos->N-1][i] = infos->s_up[3][infos->N-2][i];

        infos->s_down[0][0][i] = infos->s_down[0][1][i];
        infos->s_down[1][0][i] = infos->s_down[1][1][i];
        infos->s_down[2][0][i] = 0.0;
        infos->s_down[3][0][i] = infos->s_down[3][1][i];
    }

    //upstream
    
    for (int i = 0; i<infos->i_lead;i++){
        for (int k = 0; k<4;k++){
            infos->s_up[k][0][i] = infos->s_down[k][infos->N-1][i];
        }
    }

    //airfoil
    for(int i = infos->i_lead; i < infos->i_trail+1; i++){
        // up : j=0
        double rho = infos->s_up[0][1][i];
        double u   = infos->s_up[1][1][i] / rho;
        double v   = infos->s_up[2][1][i] / rho;
        
        // normale à la paroi (perpendiculaire à la tangente (dx, dy_updi))
        double dx = infos->dxdi[i];
        double dy = infos->dy_updi[0][i];
        double norm2 = dx*dx + dy*dy;

        double c_y0 = (u*dx + v*dy)/norm2;

        infos->s_up[0][0][i] = rho;
        infos->s_up[1][0][i] = rho *c_y0*dx ;
        infos->s_up[2][0][i] = rho *c_y0*dy ;
        infos->s_up[3][0][i] = infos->s_up[3][1][i];

        // down : j=N-1, voisin = j=N-2
        rho = infos->s_down[0][infos->N-2][i];
        u   = infos->s_down[1][infos->N-2][i] / rho;
        v   = infos->s_down[2][infos->N-2][i] / rho;
        
        dy = infos->dy_downdi[infos->N-1][i];
        norm2 = dx*dx + dy*dy;
        c_y0 = (u*dx + v*dy)/norm2;
        
        infos->s_down[0][infos->N-1][i] = rho;
        infos->s_down[1][infos->N-1][i] = rho *c_y0*dx ;
        infos->s_down[2][infos->N-1][i] = rho *c_y0*dy ;
        infos->s_down[3][infos->N-1][i] = infos->s_down[3][infos->N-2][i];
    }

    /*

    int edges[2] = {infos->i_lead, infos->i_trail};
    for (int e = 0; e < 2; e++){
        int i = edges[e];
        infos->s_up[0][0][i] = infos->s_up[0][1][i];
        infos->s_up[1][0][i] = infos->s_up[1][1][i];
        infos->s_up[2][0][i] = 0.;
        infos->s_up[3][0][i] = infos->s_up[3][1][i];

        infos->s_down[0][infos->N-1][i] = infos->s_down[0][infos->N-2][i];
        infos->s_down[1][infos->N-1][i] = infos->s_down[1][infos->N-2][i];
        infos->s_down[2][infos->N-1][i] = 0.;
        infos->s_down[3][infos->N-1][i] = infos->s_down[3][infos->N-2][i];
    }

    infos->s_up[2][0][infos->i_lead]  = 0.;
    infos->s_up[2][0][infos->i_trail] = 0.;
    infos->s_down[2][infos->N-1][infos->i_lead]  = 0.;
    infos->s_down[2][infos->N-1][infos->i_trail] = 0.;

    */
    //down stream
    for (int i = infos->i_trail+1; i<infos->M;i++){
        for (int k = 0; k<4;k++){
            infos->s_up[k][0][i] = infos->s_down[k][infos->N-1][i];
        }
    }

}

void set_bounds_s_star(struct grid_and_data *infos ){

    double rho_inf = infos->rho_inf;
    double p_inf = infos->p_inf;
    double u_inf = infos->M_inf*infos->c_sound_inf;
    double v_0 = 0;

    double rho_0, p_0, u_0, U_0;
    double rho_end, u_end, v_end, U0_end;


    for (int j = 0; j < infos->N; j++){
        // inlet conditions
        rho_0 = infos->s_star_up[0][j][1];
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        u_0 = sqrt(u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.)*(p_inf/rho_inf - p_0/rho_0));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);
        infos->s_star_up[0][j][0] = rho_0;
        infos->s_star_up[1][j][0] = rho_0*u_0;
        infos->s_star_up[2][j][0] = rho_0*v_0;
        infos->s_star_up[3][j][0] = rho_0*U_0;

        rho_0 = infos->s_star_down[0][j][1];
        p_0 = p_inf*pow(rho_0/rho_inf, infos->gamma);
        u_0 = sqrt(u_inf*u_inf + 2.*infos->gamma/(infos->gamma-1.)*(p_inf/rho_inf - p_0/rho_0));
        U_0 = p_0/((infos->gamma-1.)*rho_0) + 0.5*(u_0*u_0+v_0*v_0);
        infos->s_star_down[0][j][0] = rho_0;
        infos->s_star_down[1][j][0] = rho_0*u_0;
        infos->s_star_down[2][j][0] = rho_0*v_0;
        infos->s_star_down[3][j][0] = rho_0*U_0;
    }

    for (int j = 0; j < infos->N; j++){
        // outlet conditions
        rho_end = infos->s_star_up[0][j][infos->M-2];
        u_end = infos->s_star_up[1][j][infos->M-2]/rho_end;
        v_end = infos->s_star_up[2][j][infos->M-2]/rho_end;
        U0_end = p_inf/((infos->gamma-1.)*rho_end) + 0.5*(u_end*u_end+v_end*v_end);
        infos->s_star_up[0][j][infos->M-1] = infos->s_star_up[0][j][infos->M-2];
        infos->s_star_up[1][j][infos->M-1] = infos->s_star_up[1][j][infos->M-2];
        infos->s_star_up[2][j][infos->M-1] = infos->s_star_up[2][j][infos->M-2];
        infos->s_star_up[3][j][infos->M-1] = rho_end*U0_end;

        rho_end = infos->s_star_down[0][j][infos->M-2];
        u_end = infos->s_star_down[1][j][infos->M-2]/rho_end;
        v_end = infos->s_star_down[2][j][infos->M-2]/rho_end;
        U0_end = p_inf/((infos->gamma-1.)*rho_end) + 0.5*(u_end*u_end+v_end*v_end);
        infos->s_star_down[0][j][infos->M-1] = infos->s_star_down[0][j][infos->M-2];
        infos->s_star_down[1][j][infos->M-1] = infos->s_star_down[1][j][infos->M-2];
        infos->s_star_down[2][j][infos->M-1] = infos->s_star_down[2][j][infos->M-2];
        infos->s_star_down[3][j][infos->M-1] = rho_end*U0_end;
    }

    //  walls
    for (int i = 1; i < infos->M-1; i++){
        infos->s_star_up[0][infos->N-1][i] = infos->s_star_up[0][infos->N-2][i];
        infos->s_star_up[1][infos->N-1][i] = infos->s_star_up[1][infos->N-2][i];
        infos->s_star_up[2][infos->N-1][i] = 0.;
        infos->s_star_up[3][infos->N-1][i] = infos->s_star_up[3][infos->N-2][i];

        infos->s_star_down[0][0][i] = infos->s_star_down[0][1][i];
        infos->s_star_down[1][0][i] = infos->s_star_down[1][1][i];
        infos->s_star_down[2][0][i] = 0.;
        infos->s_star_down[3][0][i] = infos->s_star_down[3][1][i];
    }

    // frontière commune upstream
    for (int i = 0; i < infos->i_lead; i++){
        for (int k = 0; k < 4; k++)
            infos->s_star_up[k][0][i] = infos->s_star_down[k][infos->N-1][i];
    }

    //airfoil
    for(int i = infos->i_lead+1; i < infos->i_trail; i++){
        // up : j=0
        double rho = infos->s_star_up[0][1][i];
        double u   = infos->s_star_up[1][1][i] / rho;
        double v   = infos->s_star_up[2][1][i] / rho;
        
        // normale à la paroi (perpendiculaire à la tangente (dx, dy_updi))
        double dx = infos->dxdi[i];
        double dy = infos->dy_updi[0][i];
        double norm2 = dx*dx + dy*dy;

        double c_y0 = (u*dx + v*dy)/norm2;

        infos->s_star_up[0][0][i] = rho;
        infos->s_star_up[1][0][i] = rho *c_y0*dx ;
        infos->s_star_up[2][0][i] = rho *c_y0*dy ;
        infos->s_star_up[3][0][i] = infos->s_star_up[3][1][i];

        // down : j=N-1, voisin = j=N-2
        rho = infos->s_star_down[0][infos->N-2][i];
        u   = infos->s_star_down[1][infos->N-2][i] / rho;
        v   = infos->s_star_down[2][infos->N-2][i] / rho;
        
        dy = infos->dy_downdi[infos->N-1][i];
        norm2 = dx*dx + dy*dy;
        c_y0 = (u*dx + v*dy)/norm2;
        
        infos->s_star_down[0][infos->N-1][i] = rho;
        infos->s_star_down[1][infos->N-1][i] = rho *c_y0*dx ;
        infos->s_star_down[2][infos->N-1][i] = rho *c_y0*dy ;
        infos->s_star_down[3][infos->N-1][i] = infos->s_star_down[3][infos->N-2][i];
    }



    // frontière commune downstream
    for (int i = infos->i_trail+1; i < infos->M; i++){
        for (int k = 0; k < 4; k++)
            infos->s_star_up[k][0][i] = infos->s_star_down[k][infos->N-1][i];
    }
}

void init_Fs(struct grid_and_data *infos){

    int int_debug = -1;

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){

            double rho_up = infos->s_up[0][j][i];
            double u_up   = infos->s_up[1][j][i]/rho_up;
            double v_up   = infos->s_up[2][j][i]/rho_up;
            double U0_up  = infos->s_up[3][j][i]/rho_up;
            double p_up    = (infos->gamma - 1.) * rho_up * (U0_up - 0.5*(u_up*u_up + v_up*v_up));
            double H0_up = U0_up + p_up/rho_up;

            double rho_down = infos->s_down[0][j][i];
            double u_down   = infos->s_down[1][j][i]/rho_down;
            double v_down   = infos->s_down[2][j][i]/rho_down;
            double U0_down  = infos->s_down[3][j][i]/rho_down;
            double p_down   = (infos->gamma - 1.) * rho_down * (U0_down - 0.5*(u_down*u_down + v_down*v_down));
            double H0_down = U0_down + p_down/rho_down;

            infos->F_x_up[0][j][i] = infos->s_up[1][j][i];
            infos->F_x_up[1][j][i] = infos->s_up[1][j][i]*u_up + p_up;
            infos->F_x_up[2][j][i] = infos->s_up[1][j][i]*v_up;
            infos->F_x_up[3][j][i] = infos->s_up[1][j][i]*H0_up;

            infos->F_x_down[0][j][i] = infos->s_down[1][j][i];
            infos->F_x_down[1][j][i] = infos->s_down[1][j][i]*u_down + p_down;
            infos->F_x_down[2][j][i] = infos->s_down[1][j][i]*v_down;
            infos->F_x_down[3][j][i] = infos->s_down[1][j][i]*H0_down;

            infos->F_y_up[0][j][i] = infos->s_up[2][j][i];
            infos->F_y_up[1][j][i] = infos->s_up[2][j][i]*u_up;
            infos->F_y_up[2][j][i] = infos->s_up[2][j][i]*v_up + p_up;
            infos->F_y_up[3][j][i] = infos->s_up[2][j][i]*H0_up;

            infos->F_y_down[0][j][i] = infos->s_down[2][j][i];
            infos->F_y_down[1][j][i] = infos->s_down[2][j][i]*u_down;
            infos->F_y_down[2][j][i] = infos->s_down[2][j][i]*v_down + p_down;
            infos->F_y_down[3][j][i] = infos->s_down[2][j][i]*H0_down;


            // pe vérif un truc ici mais à voir avec les frontières
            for (int k = 0; k < 4;k++){
                infos->F_i_up[k][j][i] = infos->dy_updj[j][i]*infos->F_x_up[k][j][i];
                infos->F_i_down[k][j][i] = infos->dy_downdj[j][i]*infos->F_x_down[k][j][i];

                infos->F_j_up[k][j][i] = -infos->dy_updt[j][i]*infos->dxdi[i]*infos->s_up[k][j][i]- infos->dy_updi[j][i]*infos->F_x_up[k][j][i] + infos->dxdi[i]*infos->F_y_up[k][j][i];
                infos->F_j_down[k][j][i] = -infos->dy_downdt[j][i]*infos->dxdi[i]*infos->s_down[k][j][i] -infos->dy_downdi[j][i]*infos->F_x_down[k][j][i] + infos->dxdi[i]*infos->F_y_down[k][j][i];
            }
        }
    }



    for (int i = 0; i<infos->i_lead+int_debug;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_up[k][0][i] = infos->F_i_down[k][infos->N-1][i];
            infos->F_j_up[k][0][i] = infos->F_j_down[k][infos->N-1][i];            
        }
    }

    for (int i = infos->i_trail+1-int_debug; i<infos->M;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_up[k][0][i] = infos->F_i_down[k][infos->N-1][i];
            infos->F_j_up[k][0][i] = infos->F_j_down[k][infos->N-1][i];     
        }
    }
}

void init_Fs_star(struct grid_and_data *infos){

    int debug_int = -1;
    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){
            double rho_star_up = infos->s_star_up[0][j][i];
            double u_star_up   = infos->s_star_up[1][j][i]/rho_star_up;
            double v_star_up   = infos->s_star_up[2][j][i]/rho_star_up;
            double U0_star_up  = infos->s_star_up[3][j][i]/rho_star_up;
            double p_star_up   = (infos->gamma - 1.) * rho_star_up * (U0_star_up - 0.5*(u_star_up*u_star_up + v_star_up*v_star_up));
            double H0_star_up = U0_star_up + p_star_up/rho_star_up;

            double rho_star_down = infos->s_star_down[0][j][i];
            double u_star_down   = infos->s_star_down[1][j][i]/rho_star_down;
            double v_star_down   = infos->s_star_down[2][j][i]/rho_star_down;
            double U0_star_down  = infos->s_star_down[3][j][i]/rho_star_down;
            double p_star_down   = (infos->gamma - 1.) * rho_star_down * (U0_star_down - 0.5*(u_star_down*u_star_down + v_star_down*v_star_down));
            double H0_star_down = U0_star_down + p_star_down/rho_star_down;

            infos->F_x_star_up[0][j][i] = infos->s_star_up[1][j][i];
            infos->F_x_star_up[1][j][i] = infos->s_star_up[1][j][i]*u_star_up+  p_star_up;
            infos->F_x_star_up[2][j][i] = infos->s_star_up[1][j][i]*v_star_up;
            infos->F_x_star_up[3][j][i] = infos->s_star_up[1][j][i]*H0_star_up;

            infos->F_x_star_down[0][j][i] = infos->s_star_down[1][j][i];
            infos->F_x_star_down[1][j][i] = infos->s_star_down[1][j][i]*u_star_down +p_star_down;
            infos->F_x_star_down[2][j][i] = infos->s_star_down[1][j][i]*v_star_down;
            infos->F_x_star_down[3][j][i] = infos->s_star_down[1][j][i]*H0_star_down;

            infos->F_y_star_up[0][j][i] = infos->s_star_up[2][j][i];
            infos->F_y_star_up[1][j][i] = infos->s_star_up[2][j][i]*u_star_up;
            infos->F_y_star_up[2][j][i] = infos->s_star_up[2][j][i]*v_star_up +p_star_up;
            infos->F_y_star_up[3][j][i] = infos->s_star_up[2][j][i]*H0_star_up;

            infos->F_y_star_down[0][j][i] = infos->s_star_down[2][j][i];
            infos->F_y_star_down[1][j][i] = infos->s_star_down[2][j][i]*u_star_down;
            infos->F_y_star_down[2][j][i] = infos->s_star_down[2][j][i]*v_star_down+ p_star_down;
            infos->F_y_star_down[3][j][i] = infos->s_star_down[2][j][i]*H0_star_down;
            //Idem à que dans init_Fs, vérif les frontières
            for (int k = 0; k < 4;k++){
                infos->F_i_star_up[k][j][i] = infos->dy_updj[j][i]*infos->F_x_star_up[k][j][i];
                infos->F_i_star_down[k][j][i] = infos->dy_downdj[j][i]*infos->F_x_star_down[k][j][i];

                infos->F_j_star_up[k][j][i] = -infos->dy_updt[j][i]*infos->dxdi[i]*infos->s_star_up[k][j][i]-infos->dy_updi[j][i]*infos->F_x_star_up[k][j][i] + infos->dxdi[i]*infos->F_y_star_up[k][j][i];
                infos->F_j_star_down[k][j][i] = -infos->dy_downdt[j][i]*infos->dxdi[i]*infos->s_star_down[k][j][i]-infos->dy_downdi[j][i]*infos->F_x_star_down[k][j][i] + infos->dxdi[i]*infos->F_y_star_down[k][j][i];
            }
        }
    }

    for (int i = 0; i<infos->i_lead+debug_int;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_star_up[k][0][i] = infos->F_i_star_down[k][infos->N-1][i];
            infos->F_j_star_up[k][0][i] = infos->F_j_star_down[k][infos->N-1][i];            
        }
    }

    for (int i = infos->i_trail+1-debug_int; i<infos->M;i++){
        for (int k = 0; k<4;k++){
            infos->F_i_star_up[k][0][i] = infos->F_i_star_down[k][infos->N-1][i];
            infos->F_j_star_up[k][0][i] = infos->F_j_star_down[k][infos->N-1][i];    
        }
    }
}

double terme_diffusion_fct(struct grid_and_data *infos, double***s_a, double***s_b,int k,int j,int i, int j_plus_1, int j_minus_1){

    // I use s_a == s_b for every cases except the ones at the common border, where s_a = s_down and s_b = s_up = related to j+1

    double eps_min = infos->eps_min;
    double C = infos->C;
    double c = infos->c;

    int i_plus_1 = i+1;
    int i_minus_1 = i-1;

    double delta_u_i_plus_1  = s_a[1][j][i_plus_1]/ s_a[0][j][i_plus_1] - s_a[1][j][i]/ s_a[0][j][i];
    double delta_u_i_minus_1 = s_a[1][j][i]/ s_a[0][j][i] - s_a[1][j][i_minus_1]/ s_a[0][j][i_minus_1];
    double delta_u_j_plus_1  = s_b[1][j_plus_1][i]/ s_b[0][j_plus_1][i] - s_a[1][j][i]/ s_a[0][j][i];
    double delta_u_j_minus_1 = s_a[1][j][i]/ s_a[0][j][i] - s_a[1][j_minus_1][i]/ s_a[0][j_minus_1][i];

    double delta_v_i_plus_1  = s_a[2][j][i_plus_1]/ s_a[0][j][i_plus_1] - s_a[2][j][i]/ s_a[0][j][i];
    double delta_v_i_minus_1 = s_a[2][j][i]/ s_a[0][j][i] - s_a[2][j][i_minus_1]/ s_a[0][j][i_minus_1];
    double delta_v_j_plus_1  = s_b[2][j_plus_1][i]/ s_b[0][j_plus_1][i] - s_a[2][j][i]/ s_a[0][j][i];
    double delta_v_j_minus_1 = s_a[2][j][i]/ s_a[0][j][i] - s_a[2][j_minus_1][i]/ s_a[0][j_minus_1][i];

    double eps_i_plus_half  = eps_min + C*c*sqrt(delta_u_i_plus_1*delta_u_i_plus_1 + delta_v_i_plus_1*delta_v_i_plus_1);
    double eps_i_minus_half = eps_min + C*c*sqrt(delta_u_i_minus_1*delta_u_i_minus_1 + delta_v_i_minus_1*delta_v_i_minus_1);
    double eps_j_plus_half  = eps_min + C*c*sqrt(delta_u_j_plus_1*delta_u_j_plus_1 + delta_v_j_plus_1*delta_v_j_plus_1);
    double eps_j_minus_half = eps_min + C*c*sqrt(delta_u_j_minus_1*delta_u_j_minus_1 + delta_v_j_minus_1*delta_v_j_minus_1);

    double terme_diffusion_i = eps_i_plus_half*(s_a[k][j][i_plus_1]-s_a[k][j][i]) -  eps_i_minus_half*(s_a[k][j][i]-s_a[k][j][i_minus_1]);
    double terme_diffusion_j = eps_j_plus_half*(s_b[k][j_plus_1][i]-s_a[k][j][i]) -  eps_i_minus_half*(s_a[k][j][i]-s_a[k][j_minus_1][i]);

    double terme_diffusion = terme_diffusion_i + terme_diffusion_j;

    return terme_diffusion;

}

void mac_cormack_forward(struct grid_and_data *infos){

    int debug_int=-1;

    double dFidi;
    double dFjdj;
    
    double eps_min = infos->eps_min;
    double terme_diffusion;
    double terme_tot;

    //down part wall, 

    for (int i = 1; i <infos->M-1;i++ ){
        int j = 0;
        for (int k = 0; k<4;k++){
            dFidi = infos->F_i_down[k][j][i+1]-infos->F_i_down[k][j][i]; 
            dFjdj = infos->F_j_down[k][j][i] - infos->F_j_down[k][j+1][i];  // ← inversé

            // j-1 devient j+1

            //terme_diffusion = eps_min * (infos->s_down[k][j][i+1] - 2*infos->s_down[k][j][i] + infos->s_down[k][j][i-1]); // j-1 devient j+1
            terme_diffusion = terme_diffusion_fct(infos, infos->s_down,infos->s_down,k, j, i, j+1, j+1);
             
            terme_tot = -dFidi + terme_diffusion;
            infos->s_star_down[k][j][i] = infos->s_down[k][j][i]*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i] + infos->dt*terme_tot/infos->D_down_t_plus_1[j][i];        
        }
    }
    // down part
    for (int j = 1; j < infos->N-1; j++){
        for (int i = 1; i < infos->M-1;i++){
            for (int k = 0; k<4;k++){
                dFidi = infos->F_i_down[k][j][i+1]-infos->F_i_down[k][j][i];
                dFjdj = infos->F_j_down[k][j+1][i]-infos->F_j_down[k][j][i];

                //terme_diffusion = eps_min * (infos->s_down[k][j][i+1] - 2*infos->s_down[k][j][i] + infos->s_down[k][j][i-1])+ eps_min * (infos->s_down[k][j+1][i] - 2*infos->s_down[k][j][i] + infos->s_down[k][j-1][i]); 
                terme_diffusion = terme_diffusion_fct(infos, infos->s_down,infos->s_down,k, j, i, j+1, j-1);
                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_star_down[k][j][i] = infos->s_down[k][j][i]*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i] + infos->dt*terme_tot/infos->D_down_t_plus_1[j][i];
                
            }
        }
    }


    // frontière upstream  A CHANGEYYYYYYYYY
    for (int i = 1; i < infos->i_lead+debug_int;i++){
        for (int k = 0; k < 4; k++){
            int j = infos->N-1;

            dFidi = infos->F_i_down[k][j][i+1]-infos->F_i_down[k][j][i];
            dFjdj = infos->F_j_up[k][1][i]-infos->F_j_down[k][j][i];
            //dFjdj = 0.0;  // TEST
            
            terme_diffusion = terme_diffusion_fct(infos,infos->s_down, infos->s_up, k, j, i, 1, j-1 );
            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            infos->s_star_down[k][j][i] = infos->s_down[k][j][i]*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i] + infos->dt*terme_tot/infos->D_down_t_plus_1[j][i];

            infos->F_j_up[k][0][i] = infos->F_j_down[k][j][i];
            infos->F_i_up[k][0][i] = infos->F_i_down[k][j][i];
            infos->s_star_up[k][0][i] = infos->s_star_down[k][j][i];


        }
    }



    //up and downwing:

    for (int i = infos->i_lead+debug_int; i < infos->i_trail+1-debug_int;i++){
        for (int k = 0; k < 4; k++){

            int j = infos->N-1;
            //down wing, j+1 devient j-1

            dFidi = infos->F_i_down[k][j][i+1]-infos->F_i_down[k][j][i];


            //terme_diffusion = eps_min * (infos->s_down[k][j][i+1] - 2*infos->s_down[k][j][i] + infos->s_down[k][j][i-1]); 

            terme_diffusion = terme_diffusion_fct(infos, infos->s_down,infos->s_down, k, j, i, j-1, j-1);
             
            terme_tot = -dFidi  + terme_diffusion;
            infos->s_star_down[k][infos->N-1][i] = infos->s_down[k][infos->N-1][i]*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i] + infos->dt*terme_tot/infos->D_down_t_plus_1[infos->N-1][i];

            //upwing j-1 devient j+1
            j = 0;

            dFidi = infos->F_i_up[k][j][i+1]-infos->F_i_up[k][j][i];

            //terme_diffusion = eps_min * (infos->s_up[k][j][i+1] - 2*infos->s_up[k][j][i] + infos->s_up[k][j][i-1]); 
            terme_diffusion = terme_diffusion_fct(infos, infos->s_up,infos->s_up, k, j, i, j+1, j+1);
            terme_tot = -dFidi + terme_diffusion;
            //printf("forw : i = %i, terme tot down = %f", i , terme_tot);

            infos->s_star_up[k][j][i] = infos->s_up[k][j][i]*infos->D_up[j][i]/infos->D_up_t_plus_1[j][i] + infos->dt*terme_tot/infos->D_up_t_plus_1[j][i];

        }
    }

    //frontière downstream; 
    for (int i = infos->i_trail+1+debug_int; i < infos->M-1;i++){
        for (int k = 0; k < 4; k++){
            int j = infos->N-1;
            dFidi = infos->F_i_down[k][infos->N-1][i+1]-infos->F_i_down[k][infos->N-1][i];
            dFjdj = infos->F_j_up[k][1][i]-infos->F_j_down[k][infos->N-1][i];

            terme_diffusion = terme_diffusion_fct(infos,infos->s_down, infos->s_up, k, j, i, 1, j-1 );
            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            
            infos->s_star_down[k][infos->N-1][i] = infos->s_down[k][infos->N-1][i]*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i] + infos->dt*terme_tot/infos->D_down_t_plus_1[infos->N-1][i];

            infos->F_j_up[k][0][i] = infos->F_j_down[k][infos->N-1][i];
            infos->F_i_up[k][0][i] = infos->F_i_down[k][infos->N-1][i];
            infos->s_star_up[k][0][i] = infos->s_star_down[k][infos->N-1][i];
        }

    }


    //upper part exept the wall

    for (int j = 1; j < infos->N-1; j++){
        for (int i = 1; i < infos->M-1;i++){
            for (int k = 0; k<4;k++){

                dFidi = infos->F_i_up[k][j][i+1]-infos->F_i_up[k][j][i];
                dFjdj = infos->F_j_up[k][j+1][i]-infos->F_j_up[k][j][i];

                //terme_diffusion = eps_min * (infos->s_up[k][j][i+1] - 2*infos->s_up[k][j][i] + infos->s_up[k][j][i-1])+ eps_min * (infos->s_up[k][j+1][i] - 2*infos->s_up[k][j][i] + infos->s_up[k][j-1][i]);
                terme_diffusion = terme_diffusion_fct(infos, infos->s_up,infos->s_up, k, j, i, j+1, j-1);

                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_star_up[k][j][i] = infos->s_up[k][j][i]*infos->D_up[j][i]/infos->D_up_t_plus_1[j][i] + infos->dt*terme_tot/infos->D_up_t_plus_1[j][i];

            }
        }
    }

    //upper wall, j+1 devient j-1
    
    for (int i = 1; i < infos->M-1;i++){
        int j = infos->N-1;
        for (int k = 0; k<4;k++){
            dFidi = infos->F_i_up[k][j][i+1]-infos->F_i_up[k][j][i];

            //terme_diffusion = eps_min * (infos->s_up[k][j][i+1] - 2*infos->s_up[k][j][i] + infos->s_up[k][j][i-1]);
            terme_diffusion = terme_diffusion_fct(infos, infos->s_up,infos->s_up, k, j, i, j-1, j-1);

            terme_tot = -dFidi + terme_diffusion;
            infos->s_star_up[k][j][i] = infos->s_up[k][j][i]*infos->D_up[j][i]/infos->D_up_t_plus_1[j][i] + infos->dt*terme_tot/infos->D_up_t_plus_1[j][i];
        }
    }


}

void mac_cormack_backward(struct grid_and_data *infos){
    int debug_int = -1;
    double dFidi;
    double dFjdj;
    double terme_diffusion;
    double terme_tot;
    double epsilon;

    // down wall j-1 devient j+1
    for (int i = 1; i < infos->M-1;i++){
        for (int k = 0; k<4;k++){
            int j = 0;
            dFidi = infos->F_i_star_down[k][j][i]-infos->F_i_star_down[k][j][i-1];
            //dFjdj = infos->F_j_star_down[k][j+1][i] - infos->F_j_star_down[k][j][i];  


            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_down[k][j][i+1] - 2*infos->s_star_down[k][j][i] + infos->s_star_down[k][j][i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j+1, j+1);

            terme_tot = -dFidi + terme_diffusion;
            infos->s_down[k][0][i] = 0.5*(infos->s_star_down[k][0][i] + infos->s_down[k][0][i])*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i]+ 0.5*infos->dt/infos->D_down_t_plus_1[0][i] * terme_tot;
        }
    }


    //down part
    for (int j = 1; j < infos->N-1; j++){
        for (int i = 1; i < infos->M-1;i++){
            for (int k = 0; k<4;k++){

                dFidi = infos->F_i_star_down[k][j][i]-infos->F_i_star_down[k][j][i-1];
                dFjdj = infos->F_j_star_down[k][j][i]-infos->F_j_star_down[k][j-1][i];
                //dFjdj = infos->F_j_star_down[k][j+1][i] - infos->F_j_star_down[k][j][i];  // ← inversé

                epsilon = infos->eps_min;
                //terme_diffusion = epsilon * (infos->s_star_down[k][j][i+1] - 2*infos->s_star_down[k][j][i] + infos->s_star_down[k][j][i-1])+ epsilon * (infos->s_star_down[k][j-1][i] - 2*infos->s_star_down[k][j][i] + infos->s_star_down[k][j+1][i]);
                
                terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j+1, j-1);
                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_down[k][j][i] = 0.5*(infos->s_star_down[k][j][i] + infos->s_down[k][j][i])*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i]+ 0.5*infos->dt/infos->D_down_t_plus_1[j][i] * terme_tot;

            }
        }
    }

    // frontière upstream 
    for (int i = 1; i < infos->i_lead+debug_int;i++){
        int j = infos->N-1;
        for (int k = 0; k < 4; k++){
            dFidi = infos->F_i_star_down[k][j][i]-infos->F_i_star_down[k][j][i-1];
            dFjdj = infos->F_j_star_down[k][j][i]-infos->F_j_star_up[k][1][i];
            //dFjdj = 0.0;  // TEST

            epsilon = infos->eps_min;
            terme_diffusion = epsilon * (infos->s_star_down[k][j][i+1] - 2*infos->s_star_down[k][j][i] + infos->s_star_down[k][j][i-1])+ epsilon * (infos->s_star_down[k][j-1][i] - 2*infos->s_star_down[k][j][i] + infos->s_star_up[k][1][i]);

            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            infos->s_down[k][j][i] = 0.5*(infos->s_star_down[k][j][i] + infos->s_down[k][j][i])*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i] + 0.5*infos->dt/infos->D_down_t_plus_1[j][i] * terme_tot;

            infos->F_j_star_up[k][0][i] = infos->F_j_star_down[k][j][i];
            infos->F_i_star_up[k][0][i] = infos->F_i_star_down[k][j][i];
            //infos->s_up[k][0][i] = infos->s_down[k][infos->N-1][i];
        }
    }

    //up and down wing:
    
    for (int i = infos->i_lead+debug_int; i < infos->i_trail+1-debug_int;i++){
        for (int k = 0; k < 4; k++){
            int j = infos->N-1;
            //down wing, j+1 devient j-1

            dFidi = infos->F_i_star_down[k][j][i]-infos->F_i_star_down[k][j][i-1];



            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_down[k][j][i+1] - 2*infos->s_star_down[k][j][i] + infos->s_star_down[k][j][i-1]);
            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_down,infos->s_star_down, k, j, i, j-1, j-1);

            terme_tot = -dFidi + terme_diffusion;
            infos->s_down[k][j][i] = 0.5*(infos->s_star_down[k][j][i] + infos->s_down[k][j][i])*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i] + 0.5*infos->dt/infos->D_down_t_plus_1[j][i] * terme_tot;

            //upwing j-1 devient j+1
            j = 0;
            dFidi = infos->F_i_star_up[k][j][i]-infos->F_i_star_up[k][j][i-1];

            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_up[k][j][i+1] - 2*infos->s_star_up[k][j][i] + infos->s_star_up[k][j][i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,  infos->s_star_up, k, j, i, j+1, j+1);
            terme_tot = -dFidi + terme_diffusion;

            infos->s_up[k][j][i] = 0.5*(infos->s_star_up[k][j][i] + infos->s_up[k][j][i])*infos->D_up[j][i]/infos->D_up_t_plus_1[j][i]+ 0.5*infos->dt/infos->D_up_t_plus_1[j][i] * terme_tot;

        }

    }
    //frontière downstream;
    for (int i = infos->i_trail+1-debug_int; i < infos->M-1;i++){
        for (int k = 0; k < 4; k++){
            int j = infos->N-1;
            dFidi = infos->F_i_star_down[k][j][i]-infos->F_i_star_down[k][j][i-1];
            dFjdj = infos->F_j_star_down[k][j][i]-infos->F_j_star_up[k][1][i];
            //dFjdj = 0.0;  // TEST

            epsilon = infos->eps_min;
            terme_diffusion = epsilon * (infos->s_star_down[k][j][i+1] - 2*infos->s_star_down[k][j][i] + infos->s_star_down[k][j][i-1])+ epsilon * (infos->s_star_down[k][j-1][i] - 2*infos->s_star_down[k][j][i] + infos->s_star_up[k][1][i]);

            terme_tot = -(dFidi + dFjdj) + terme_diffusion;
            infos->s_down[k][j][i] = 0.5*(infos->s_star_down[k][j][i] + infos->s_down[k][j][i])*infos->D_down[j][i]/infos->D_down_t_plus_1[j][i] + 0.5*infos->dt/infos->D_down_t_plus_1[j][i] * terme_tot;

            infos->F_j_star_up[k][0][i] = infos->F_j_star_down[k][j][i];
            infos->F_i_star_up[k][0][i] = infos->F_i_star_down[k][j][i];
            //infos->s_up[k][0][i] = infos->s_down[k][j][i];

        }
    }
    //upper part
    for (int j = 1; j < infos->N-1; j++){
        for (int i = 1; i < infos->M-1;i++){
            for (int k = 0; k<4;k++){

                dFidi = infos->F_i_star_up[k][j][i]-infos->F_i_star_up[k][j][i-1];
                dFjdj = infos->F_j_star_up[k][j][i]-infos->F_j_star_up[k][j-1][i];

                epsilon = infos->eps_min;
                //terme_diffusion = epsilon * (infos->s_star_up[k][j][i+1] - 2*infos->s_star_up[k][j][i] + infos->s_star_up[k][j][i-1])+ epsilon * (infos->s_star_up[k][j+1][i] - 2*infos->s_star_up[k][j][i] + infos->s_star_up[k][j-1][i]);
                terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,infos->s_star_up, k, j, i, j+1, j-1);
                terme_tot = -(dFidi + dFjdj) + terme_diffusion;
                infos->s_up[k][j][i] = 0.5*(infos->s_star_up[k][j][i] + infos->s_up[k][j][i])*infos->D_up[j][i]/infos->D_up_t_plus_1[j][i] + 0.5*infos->dt/infos->D_up_t_plus_1[j][i] * terme_tot;
            }
        }

    }

    // upper wall j+1 devient j-1

    for (int i = 1; i < infos->M-1;i++){
        int j = infos->N-1;
        for (int k = 0; k<4;k++){

            dFidi = infos->F_i_star_up[k][j][i]-infos->F_i_star_up[k][j][i-1];

            epsilon = infos->eps_min;
            //terme_diffusion = epsilon * (infos->s_star_up[k][j][i+1] - 2*infos->s_star_up[k][j][i] + infos->s_star_up[k][j][i-1]);

            terme_diffusion = terme_diffusion_fct(infos, infos->s_star_up,infos->s_star_up, k, j, i, j-1, j+1);
        
            terme_tot = -dFidi  + terme_diffusion;
            infos->s_up[k][j][i] = 0.5*(infos->s_star_up[k][j][i] + infos->s_up[k][j][i])*infos->D_up[j][i]/infos->D_up_t_plus_1[j][i] + 0.5*infos->dt/infos->D_up_t_plus_1[j][i] * terme_tot;
        }
    }
}

void solver_one_time_step(struct grid_and_data *infos){

    init_Fs(infos);
    mac_cormack_forward(infos);
    set_bounds_s_star(infos); 
    init_Fs_star(infos);
    mac_cormack_backward(infos);
    set_bounds_s(infos);

}

void solver_all_steps(struct grid_and_data *infos){

    set_x_grid(infos);
    leading_and_trailing_edge_finder(infos);   
    init_y_grid(infos);
    update_y_grid(infos);
    set_init_s(infos);
    set_bounds_s(infos);

    FILE * file = generate_file(infos);
    write_first_two_lines(infos, file); 

    int steps = 0;


    while(infos->t < infos->t_end){
        solver_one_time_step(infos);
        if (steps == infos->step_max){
            write_one_time_step(infos, file); 
            steps = 0;
        }
        update_y_grid(infos);
        infos->t += infos->dt;

        steps +=1;
    }

    fprintf(file, "end \n");
    fclose(file);

}

#endif 