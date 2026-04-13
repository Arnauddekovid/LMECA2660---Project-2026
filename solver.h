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



void set_init_s(struct grid_and_data *infos){

    double rho_inf = infos->rho_inf;
    double u_inf = infos->M_inf*infos->c_sound_inf;
    double v_inf = 0.;
    double U_0 = infos->p_inf/((infos->gamma-1)*rho_inf) + 0.5*(u_inf*u_inf+v_inf*v_inf);
    
    for (int j =0; j <infos->N;j++){
        for (int i = 0 ;  i < infos->M;i++){
            
            infos->s_up[0][j][i] = rho_inf;
            infos->s_up[1][j][i] = rho_inf*u_inf;
            infos->s_up[2][j][i] = rho_inf*v_inf;
            infos->s_up[3][j][i] = rho_inf*U_0;            
            
            infos->s_down[0][j][i] = rho_inf;
            infos->s_down[1][j][i] = rho_inf*u_inf;
            infos->s_down[2][j][i] = rho_inf*v_inf;
            infos->s_down[3][j][i] = rho_inf*U_0;
        }
    }


}

void set_bounds_s(struct grid_and_data *infos ){

    double rho_inf = infos->rho_inf;
    double p_inf = infos->p_inf;
    double u_inf = infos->M_inf*infos->c_sound_inf;
    double v_inf = 0.;
    double U_0 = p_inf/((infos->gamma-1.)*rho_inf) + 0.5*(u_inf*u_inf+v_inf*v_inf);

    double rho_end;
    double u_end;
    double v_end;
    double U0_end;

    for (int j = 0; j< infos->N; j++){
        //inlet conditions
        infos->s_up[0][j][0] = rho_inf;
        infos->s_up[1][j][0] = rho_inf*u_inf;
        infos->s_up[2][j][0] = rho_inf*v_inf;
        infos->s_up[3][j][0] = rho_inf*U_0;

        infos->s_down[0][j][0]= rho_inf;
        infos->s_down[1][j][0]= rho_inf*u_inf;
        infos->s_down[2][j][0]= rho_inf*v_inf;
        infos->s_down[3][j][0]= rho_inf*U_0;

        //outlet conditions

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

    for (int i = 0; i <infos->M; i++){

        // le wall est en j-1 en up
        rho_wall = infos->s_up[0][infos->N-2][i];
        u_wall = infos->s_up[1][infos->N-2][i]/rho_wall;
        v_wall = infos->s_up[2][infos->N-2][i]/rho_wall;
        U_wall = p_inf/((infos->gamma -1)*rho_wall) + 0.5*(u_wall*u_wall+v_wall*v_wall);

        infos->s_up[0][infos->N-1][i] = rho_wall;
        infos->s_up[1][infos->N-1][i] = rho_wall*u_wall;
        infos->s_up[2][infos->N-1][i] = rho_wall*v_wall;
        infos->s_up[3][infos->N-1][i] = rho_wall*U_wall;

        //le wall est en j = 0 en down
        rho_wall = infos->s_down[0][1][i];
        u_wall = infos->s_down[1][1][i]/rho_wall;
        v_wall = infos->s_up[2][1][i]/rho_wall;
        U_wall = p_inf/((infos->gamma -1)*rho_wall) + 0.5*(u_wall*u_wall+v_wall*v_wall);

        infos->s_down[0][0][i] = rho_wall;
        infos->s_down[1][0][i] = rho_wall*u_wall;
        infos->s_down[2][0][i] = rho_wall*v_wall;
        infos->s_down[3][0][i] = rho_wall*U_wall;


    }
}

void set_bounds_s_star(struct grid_and_data *infos){


    double rho_inf = infos->rho_inf;
    double p_inf = infos->p_inf;
    
    double u_inf = infos->M_inf*infos->c_sound_inf;
    double v_inf = 0.;
    double U_0 = p_inf/((infos->gamma-1.)*rho_inf) + 0.5*(u_inf*u_inf + v_inf*v_inf);

    double rho_end;
    double u_end;
    double v_end;
    double U0_end;

    for (int j = 0; j< infos->N; j++){
        //inlet conditions
        infos->s_star_up[0][j][0] = rho_inf;
        infos->s_star_up[1][j][0] = rho_inf*u_inf;
        infos->s_star_up[2][j][0] = rho_inf*v_inf;
        infos->s_star_up[3][j][0] = rho_inf*U_0;

        infos->s_star_down[0][j][0]= rho_inf;
        infos->s_star_down[1][j][0]= rho_inf*u_inf;
        infos->s_star_down[2][j][0]= rho_inf*v_inf;
        infos->s_star_down[3][j][0]= rho_inf*U_0;
        //outlet conditions

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

    // wall dumby bumd conditions ( on prend pas encore en compte les trucs courbés)
    // à priori c'es u dx - v dy = 0 mais à vérifier

    double rho_wall; 
    double u_wall;
    double v_wall;
    double U_wall;

    for (int i = 0; i <infos->M; i++){

        // le wall est en j-1 en up
        rho_wall = infos->s_star_up[0][infos->N-2][i];
        u_wall = infos->s_star_up[1][infos->N-2][i]/rho_wall;
        v_wall = infos->s_star_up[2][infos->N-2][i]/rho_wall;
        U_wall = p_inf/((infos->gamma -1)*rho_wall) + 0.5*(u_wall*u_wall+v_wall*v_wall);

        infos->s_star_up[0][infos->N-1][i] = rho_wall;
        infos->s_star_up[1][infos->N-1][i] = rho_wall*u_wall;
        infos->s_star_up[2][infos->N-1][i] = rho_wall*v_wall;
        infos->s_star_up[3][infos->N-1][i] = rho_wall*U_wall;

        //le wall est en j = 0 en down
        rho_wall = infos->s_star_down[0][1][i];
        u_wall = infos->s_star_down[1][1][i]/rho_wall;
        v_wall = infos->s_star_down[2][1][i]/rho_wall;
        U_wall = p_inf/((infos->gamma -1)*rho_wall) + 0.5*(u_wall*u_wall+v_wall*v_wall);

        infos->s_star_down[0][0][i] = rho_wall;
        infos->s_star_down[1][0][i] = rho_wall*u_wall;
        infos->s_star_down[2][0][i] = rho_wall*v_wall;
        infos->s_star_down[3][0][i] = rho_wall*U_wall;


    }
}


void solver_one_time_step(struct grid_and_data *infos){
    // initialisation of the F_x and F_y stuffs
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
            infos->F_x_up[1][j][i] = infos->s_up[1][j][i]*u_up+  p_up;
            infos->F_x_up[2][j][i] = infos->s_up[1][j][i]*v_up;
            infos->F_x_up[3][j][i] = infos->s_up[1][j][i]*H0_up;

            infos->F_x_down[0][j][i] = infos->s_down[1][j][i];
            infos->F_x_down[1][j][i] = infos->s_down[1][j][i]*u_down +p_down;
            infos->F_x_down[2][j][i] = infos->s_down[1][j][i]*v_down;
            infos->F_x_down[3][j][i] = infos->s_down[1][j][i]*H0_down;

            infos->F_y_up[0][j][i] = infos->s_up[2][j][i];
            infos->F_y_up[1][j][i] = infos->s_up[2][j][i]*u_up;
            infos->F_y_up[2][j][i] = infos->s_up[2][j][i]*v_up +p_up;
            infos->F_y_up[3][j][i] = infos->s_up[2][j][i]*H0_up;

            infos->F_y_down[0][j][i] = infos->s_down[2][j][i];
            infos->F_y_down[1][j][i] = infos->s_down[2][j][i]*u_down;
            infos->F_y_down[2][j][i] = infos->s_down[2][j][i]*v_down+ p_down;
            infos->F_y_down[3][j][i] = infos->s_down[2][j][i]*H0_down;

            //pas vérif
            for (int k = 0; k < 4;k++){
                infos->F_i_up[k][j][i] = infos->dy_updj[j][i]*infos->F_x_up[k][j][i];
                infos->F_i_down[k][j][i] = infos->dy_downdj[j][i]*infos->F_x_down[k][j][i];

                infos->F_j_up[k][j][i] = -infos->dy_updt[j][i]*infos->dxdi[i]*infos->s_up[k][j][i]- infos->dy_updi[j][i]*infos->F_x_up[k][j][i] + infos->dxdi[i]*infos->F_y_up[k][j][i];
                infos->F_j_down[k][j][i] = -infos->dy_downdt[j][i]*infos->dxdi[i]*infos->s_down[k][j][i] -infos->dy_downdi[j][i]*infos->F_x_down[k][j][i] + infos->dxdi[i]*infos->F_y_down[k][j][i];
            }
        }
    }

    for (int k = 0; k < 4; k++){
        for (int j = 0; j < infos->N; j++){
            for (int i = 0; i < infos->M; i++){
                infos->s_star_up[k][j][i]   = infos->s_up[k][j][i];
                infos->s_star_down[k][j][i] = infos->s_down[k][j][i];
            }
        }
    }

    //forward
    for (int j = 0; j < infos->N-1; j++){
        for (int i =0; i < infos->M-1;i++){
            for (int k = 0; k<4;k++){
                double dFidi_up = infos->F_i_up[k][j][i+1]-infos->F_i_up[k][j][i];
                double dFjdj_up = infos->F_j_up[k][j+1][i]-infos->F_j_up[k][j][i];
                double terme_diffusion_up = 0.; /// A CHANGEYYY
                double terme_up = -(dFidi_up + dFjdj_up) + terme_diffusion_up;

                double dFidi_down = infos->F_i_down[k][j][i+1]-infos->F_i_down[k][j][i];
                double dFjdj_down = infos->F_j_down[k][j+1][i]-infos->F_j_down[k][j][i];
                double terme_diffusion_down = 0.; /// A CHANGEYYY
                double terme_down = -(dFidi_down + dFjdj_down) + terme_diffusion_down;

                infos->s_star_up[k][j][i] = infos->s_up[k][j][i] - infos->dt*terme_up/infos->D_up[j][i];
                infos->s_star_down[k][j][i] = infos->s_down[k][j][i] - infos->dt*terme_down/infos->D_down[j][i];

            }
        }
    }

    set_bounds_s_star(infos);

    // corrector part

    for (int j = 0; j < infos->N; j++){
        for (int i = 0; i < infos->M; i++){

            double rho_star_up = infos->s_star_up[0][j][i];
            double u_star_up   = infos->s_star_up[1][j][i]/rho_star_up;
            double v_star_up   = infos->s_star_up[2][j][i]/rho_star_up;
            double U0_star_up  = infos->s_star_up[3][j][i]/rho_star_up;
            double p_star_up    = (infos->gamma - 1.) * rho_star_up * (U0_star_up - 0.5*(u_star_up*u_star_up + v_star_up*v_star_up));
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
            //pas vérif
            for (int k = 0; k < 4;k++){
                infos->F_i_star_up[k][j][i] = infos->dy_updj[j][i]*infos->F_x_star_up[k][j][i];
                infos->F_i_star_down[k][j][i] = infos->dy_downdj[j][i]*infos->F_x_star_down[k][j][i];

                infos->F_j_star_up[k][j][i] = -infos->dy_updt[j][i]*infos->dxdi[i]*infos->s_star_up[k][j][i]-infos->dy_updi[j][i]*infos->F_x_star_up[k][j][i] + infos->dxdi[i]*infos->F_y_star_up[k][j][i];
                infos->F_j_star_down[k][j][i] = -infos->dy_downdt[j][i]*infos->dxdi[i]*infos->s_star_down[k][j][i]-infos->dy_downdi[j][i]*infos->F_x_star_down[k][j][i] + infos->dxdi[i]*infos->F_y_star_down[k][j][i];

            }
        }
    }

    for (int j = 1; j < infos->N; j++){
        for (int i = 1; i < infos->M;i++){
            for (int k = 0; k<4;k++){

                double dFidi_up = infos->F_i_star_up[k][j][i]-infos->F_i_star_up[k][j][i-1];
                double dFjdj_up = infos->F_j_star_up[k][j][i]-infos->F_j_star_up[k][j-1][i];
                double terme_diffusion_up = 0.; /// A CHANGEYYY
                double terme_up = -(dFidi_up + dFjdj_up) + terme_diffusion_up;

                double dFidi_down = infos->F_i_star_down[k][j][i]-infos->F_i_star_down[k][j][i-1];
                double dFjdj_down = infos->F_j_star_down[k][j][i]-infos->F_j_star_down[k][j-1][i];
                double terme_diffusion_down = 0.; /// A CHANGEYYY
                double terme_down = -(dFidi_down + dFjdj_down) + terme_diffusion_down;

                infos->s_up[k][j][i] = 0.5*(infos->s_star_up[k][j][i] + infos->s_up[k][j][i])+ 0.5*infos->dt/infos->D_up[j][i] * terme_up;
                infos->s_down[k][j][i] = 0.5*(infos->s_star_down[k][j][i] + infos->s_down[k][j][i]) + 0.5*infos->dt/infos->D_down[j][i] * terme_up;

            }
        }
    }

    set_bounds_s(infos);

    infos->t +=infos->dt;

}

void solver_all_steps(struct grid_and_data *infos){

    set_init_s(infos);
    FILE * file = generate_file(infos);
    write_first_line(infos, file); // writes the x infos

    while(infos->t < infos->t_end){
        solver_one_time_step(infos);
        write_one_time_step(infos, file);
    }
    fprintf(file, "end \n");
    fclose(file);
}
#endif 