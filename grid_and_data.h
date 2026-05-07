//
//  grid_and_data.h
//  MECA2660 projet 2026

#ifndef grid_and_data_h
#define grid_and_data_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <stdbool.h>  

struct grid_and_data{

    //pysical data 
    double M_inf_init;
    double M_inf;
    double grow_rate_M_inf;
    double gamma;
    double R_fluid;
    double c_sound_inf;
    double T_inf;
    double T_0_inf;
    double p_inf;
    double p_0_inf;
    double rho_inf;

    double C; // le petit facteur pr la diffusion
    double eps_min;

    //geometric data
    double c;
    double R; // the radius of the circle
    double H;
    double y_c;
    double x_c;

    double L;
    double L_u;


    // points i formations
    int M;
    int i_1;
    int i_2;
    int i_lead;
    int i_trail;
    int N;

    //iterations
    double t;
    double dt;
    double t_end;
    double step_max;

    // some parameter we must choose carefuly
    double beta_1;
    double beta_3;
    double beta_minus;
    double beta_plus;
    double beta;

    double AOA;
    double mu;
    double omega;

    // grids
    double* x;  
    double**y_up;
    double**y_up_t_plus_1; // the grid at the t + dt in order to get the tempo dervative
    double**y_down;
    double**y_down_t_plus_1; 




    //derivatives and other stuff

    double** D_up;
    double** D_down;
    double** D_up_t_plus_1; // à changer quand ça bougera au cours du temps
    double** D_down_t_plus_1; // à changer quand ça bougera au cours du temps

    double* dxdi;
    double** dy_updi;
    double** dy_updj;
    double** dy_updt;
    double** dy_downdi;
    double** dy_downdj;
    double** dy_downdt;
    //no dxdj and no dxdt


    // physical dimensions in terms of s and F
    
    double ***s_up;
    double ***s_down;
    double ***s_star_up;
    double ***s_star_down;

    double ***F_x_up;
    double ***F_x_down;
    double ***F_x_star_up; 
    double ***F_x_star_down; 

    double ***F_y_up;
    double ***F_y_down;
    double ***F_y_star_up; 
    double ***F_y_star_down; 

    double ***F_i_up;
    double ***F_i_down;
    double ***F_i_star_up; 
    double ***F_i_star_down; 

    double ***F_j_up;
    double ***F_j_down;
    double ***F_j_star_up; 
    double ***F_j_star_down; 
};

void grow_L_to_7(struct grid_and_data *infos){
    infos->L = 7.;
    infos->beta_3 = 1.93668189;
}

void leading_and_trailing_edge_finder(struct grid_and_data *infos){
    int i = 0;

    while (infos->x[i] < 0.){
        
        i +=1;
    }
    infos->i_lead = i-1;
    while (infos->x[i] <1.){
        
        i +=1;
    }
    infos->i_trail = i;

}

double f_plus(struct grid_and_data*infos, int index){   
    
    double x_i = infos->x[index];
    double value_space = sqrt(-(x_i-infos->x_c)*(x_i-infos->x_c)+infos->R*infos->R)- sqrt(-infos->x_c*infos->x_c+infos->R*infos->R);
    double value_time = (1.-2.*x_i/infos->c)*infos->mu*infos->c*sin(infos->omega*infos->t);
    return (value_space +value_time)  /2.;
}

double f_minus(struct grid_and_data*infos, int index){
    double x_i = infos->x[index];
    double value_space = -sqrt(-(x_i-infos->x_c)*(x_i-infos->x_c)+infos->R*infos->R)+ sqrt(-infos->x_c*infos->x_c+infos->R*infos->R);
    double value_time  = (1.-2.*x_i/infos->c)*infos->mu*infos->c*sin(infos->omega*infos->t);
    return ( value_space+value_time) /2.;
}

double fonction_time_1(struct grid_and_data *infos, double time){

    return infos->mu*infos->c*sin(infos->omega*time) ;
}

void set_x_grid(struct grid_and_data *infos){

    // partie pour set l'axe x et dxdi
    // upstream
    for (int i= 0; i<infos->i_1; i++){
        infos->x[i] = (infos->L_u - infos->c/2)* tanh(infos->beta_1*i/infos->i_1)/tanh(infos->beta_1) - infos->L_u;

        double arg = (double)infos->beta_1/infos->i_1;

        double numerator = (infos->L_u - infos->c/2)*arg;
        double denominator = tanh(infos->beta_1)*cosh(arg*i)*cosh(arg*i);
        infos->dxdi[i] = numerator/denominator;
    }

    //midstream
    for (int i = infos->i_1; i< infos->i_2; i++){
        infos->x[i] = ( (infos->L_u - infos->c/2) + 2*infos->c*(i-infos->i_1)/(infos->i_2 - infos->i_1))- infos->L_u;
        infos->dxdi[i] = 2*infos->c/(infos->i_2 - infos->i_1);
    }
    // downstream
    double terme1 = (infos->L_u + 3*infos->c/2)/infos->L;

    for (int i = infos->i_2; i<infos->M; i ++){

        double term1 = (infos->L_u + 3*infos->c/2)/infos->L;
        double term2 = 1. - term1;

        int truc_num = (i - infos->i_2);
        int truc_den = (infos->M-1-infos->i_2);
        double truc = (double)truc_num / (double)truc_den;
        double f_i = 1.-truc;

        double arg = infos->beta_3* (1.- (double)truc);

        double term3 =  1. - ( tanh ( arg)/tanh(infos->beta_3));

        infos->x[i] = infos->L*(term1 + term2*term3)-infos->L_u;

        double fprim_i = -(double)1./(infos->M-1.-infos->i_2);

        double numerator = -infos->L*(term2)*infos->beta_3*fprim_i;
        double denominator = tanh(infos->beta_3)*cosh(infos->beta_3*f_i)*cosh(infos->beta_3*f_i);

        infos->dxdi[i] = numerator/denominator;
        
    }

}

void init_y_grid(struct grid_and_data *infos){

    //on calcule les points en y, on calcule les dérivées dans update_y_grid
    // trouve les indexes tq x = 0 et x=1
    /*
    double time_coef = fonction_time_1(infos, infos->t ) ;
    double arg_up;
    double arg_down;

    double coef_up;
    double coef_down;

    for (int j =0; j <infos->N;j++){
        arg_up =  1. - (double)j/(infos->N-1);
        arg_down = (double)j/(infos->N-1);

        coef_up =   1. - tanh(infos->beta_plus*arg_up)/tanh(infos->beta);
        coef_down = 1. - tanh(infos->beta_minus* arg_down)/tanh(infos->beta);
        
        for (int i = 0; i < infos->i_lead; i++){
            infos->y_up_t_plus_1[j][i]   = (infos->H/2- time_coef)*coef_up + time_coef;
            infos->y_down_t_plus_1[j][i] =  -(infos->H/2 +time_coef)*coef_down +time_coef;
        }

        for (int i = infos->i_lead; i < infos->i_trail+1; i++){
            infos->y_up_t_plus_1[j][i] = infos->H/2*(coef_up + arg_up*f_plus(infos, i));
            infos->y_down_t_plus_1[j][i] = -infos->H/2*(coef_down - arg_down*f_minus(infos, i));
        }

        for(int i = infos->i_trail+1;i<infos->M;i++){
            infos->y_up_t_plus_1[j][i] =  (infos->H/2 +time_coef)*coef_up  - time_coef;
            infos->y_down_t_plus_1[j][i] = -(infos->H/2-time_coef)*coef_down -time_coef;
        }
    }
    */

    int i_lead  = infos->i_lead;
    int i_trail = infos->i_trail;
    int di = 8;

    double time_coef          = fonction_time_1(infos, infos->t);
    double time_coef_t_plus_1 = fonction_time_1(infos, infos->t + infos->dt);

    for (int j = 0; j < infos->N; j++){

        double arg_up   = 1. - (double)j / (infos->N-1); 
        double arg_down = (double)j / (infos->N-1);       

        double coef_up   = 1. - tanh(infos->beta_plus *arg_up)   / tanh(infos->beta);
        double coef_down = tanh(infos->beta_minus*arg_down)  / tanh(infos->beta) -1.;

        for (int i = 0; i < infos->M; i++){
            infos->y_up[j][i]   = infos->y_up_t_plus_1[j][i];
            infos->y_down[j][i] = infos->y_down_t_plus_1[j][i];
        }

        // upstream 
        for (int i = 0; i < i_lead; i++){
            infos->y_up_t_plus_1[j][i]   = (infos->H/2 - time_coef_t_plus_1)*coef_up   + time_coef_t_plus_1;
            infos->y_down_t_plus_1[j][i] = (infos->H/2 + time_coef_t_plus_1)*coef_down + time_coef_t_plus_1;
        }

        // leading edge smoothin
        for (int i = i_lead; i < i_lead + di; i++){
            double step  = (double)(i - i_lead) / (double)di;
            double scale = step*step*(3. - 2.*step);  // 0 → 1

            // Valeur upstream (ce qu'on avait juste avant i_lead)
            double y_up_upstream   = (infos->H/2 - time_coef)*coef_up   + time_coef;
            double y_down_upstream = (infos->H/2 + time_coef)*coef_down + time_coef;

            // Valeur airfoil complet (ce qu'on aura à i_lead+di)
            double y_up_airfoil   = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *f_plus(infos, i);
            double y_down_airfoil = (infos->H/2)*coef_down + (infos->H/2)*arg_down*f_minus(infos, i);

            // Interpolation cubique entre les deux
            infos->y_up_t_plus_1[j][i]   = (1.-scale)*y_up_upstream   + scale*y_up_airfoil;
            infos->y_down_t_plus_1[j][i] = (1.-scale)*y_down_upstream + scale*y_down_airfoil;
        }

        // arifoil central
        for (int i = i_lead + di; i <= i_trail - di; i++){
            double df_plusdi  = -0.5*(infos->x[i]-infos->x_c)/ sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c));
            double df_minusdi = -df_plusdi;

            infos->y_up_t_plus_1[j][i]   = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *f_plus(infos, i);
            infos->y_down_t_plus_1[j][i] = (infos->H/2)*coef_down + (infos->H/2)*arg_down*f_minus(infos, i);;


        }

        // Trailing edge smoothing
        for (int i = i_trail - di; i <= i_trail; i++){
            double step  = (double)(i_trail - i) / (double)di;
            double scale = step*step*(3. - 2.*step);  // 1 → 0

            // Valeur downstream (ce qu'on aura juste après i_trail)
            double y_up_downstream   = (infos->H/2 + time_coef)*coef_up   - time_coef;
            double y_down_downstream = (infos->H/2 - time_coef)*coef_down - time_coef;

            // Valeur airfoil complet (ce qu'on avait à i_trail-di)
            double y_up_airfoil   = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *f_plus(infos, i);
            double y_down_airfoil = (infos->H/2)*coef_down + (infos->H/2)*arg_down*f_minus(infos, i);

            // Interpolation cubique entre airfoil et downstream
            infos->y_up_t_plus_1[j][i]   = scale*y_up_airfoil   + (1.-scale)*y_up_downstream;
            infos->y_down_t_plus_1[j][i] = scale*y_down_airfoil + (1.-scale)*y_down_downstream;
        }


        // downstream
        for (int i = i_trail+1; i < infos->M; i++){
            infos->y_up_t_plus_1[j][i]   =  (infos->H/2 + time_coef_t_plus_1)*coef_up   - time_coef_t_plus_1;
            infos->y_down_t_plus_1[j][i] =  (infos->H/2 - time_coef_t_plus_1)*coef_down - time_coef_t_plus_1;

        }
    }
}

void update_y_grid_assymetrie(struct grid_and_data *infos){

    infos->M_inf = infos->M_inf_init*(1. + infos->grow_rate_M_inf*infos->t);

    int i_lead  = infos->i_lead;
    int i_trail = infos->i_trail;
    int di = 1;

    double time_coef        =  fonction_time_1(infos, infos->t);
    double time_coef_t_plus_1 = fonction_time_1(infos, infos->t + infos->dt);

    double coef_coef = 2.;

    for (int j = 0; j < infos->N; j++){

        double arg_up   = 1. - (double)j / (infos->N-1); 
        double arg_down = (double)j / (infos->N-1);       

        double coef_up   = 1. - tanh(infos->beta_plus *arg_up)   / tanh(infos->beta);
        double coef_down = tanh(infos->beta_minus*arg_down)  / tanh(infos->beta) -1.;

        double h_j_up   = coef_up;   
        double h_j_down = coef_down;

        double numerator_up   = infos->beta_plus  / (double)(infos->N-1);
        double denominator_up = tanh(infos->beta)*cosh(infos->beta_plus*arg_up)*cosh(infos->beta_plus*arg_up);
        double dh_updj_coef   = (infos->H/2)*numerator_up / denominator_up;

        double numerator_down   = infos->beta_minus / (double)(infos->N-1);
        double denominator_down = tanh(infos->beta)*cosh(infos->beta_minus*arg_down)*cosh(infos->beta_minus*arg_down);
        double dh_downdj_coef   = (infos->H/2)*numerator_down / denominator_down;


        for (int i = 0; i < infos->M; i++){
            infos->y_up[j][i]   = infos->y_up_t_plus_1[j][i];
            infos->y_down[j][i] = infos->y_down_t_plus_1[j][i];
            infos->dy_updi[j][i]   = 0.;
            infos->dy_downdi[j][i] = 0.;
        }

        // upstream 
        for (int i = 0; i < i_lead; i++){
            infos->y_up_t_plus_1[j][i]   = (infos->H/2 - coef_coef*time_coef_t_plus_1)*coef_up   + coef_coef*time_coef_t_plus_1;
            infos->y_down_t_plus_1[j][i] = (infos->H/2 + coef_coef*time_coef_t_plus_1)*coef_down + coef_coef*time_coef_t_plus_1;

            infos->dy_updj[j][i]   = (infos->H/2 - time_coef)*numerator_up   / denominator_up;
            infos->dy_downdj[j][i] = (infos->H/2 + time_coef)*numerator_down / denominator_down;
        }

        // Airfoil central
        for (int i = i_lead ; i <= i_trail; i++){
            double y_c_here          = time_coef          * (1. - 2.*infos->x[i] / infos->c);
            double y_c_here_t_plus_1 = time_coef_t_plus_1 * (1. - 2.*infos->x[i] / infos->c);
            double dy_c_here_di      = -2.*time_coef / infos->c * infos->dxdi[i];

            double df_plusdi  = -0.5*(infos->x[i]-infos->x_c) / sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c));
            double df_minusdi = -df_plusdi;

            infos->y_up_t_plus_1[j][i]   = (infos->H/2 - y_c_here_t_plus_1)*coef_up   + y_c_here_t_plus_1 + (infos->H/2)*arg_up  *f_plus(infos, i);
            infos->y_down_t_plus_1[j][i] = (infos->H/2 + y_c_here_t_plus_1)*coef_down + y_c_here_t_plus_1 + (infos->H/2)*arg_down*f_minus(infos, i);

            infos->dy_updj[j][i]   = (infos->H/2 - y_c_here)*numerator_up   / denominator_up   + dh_updj_coef  *f_plus(infos, i);
            infos->dy_downdj[j][i] = (infos->H/2 + y_c_here)*numerator_down / denominator_down - dh_downdj_coef*f_minus(infos, i);

            infos->dy_updi[j][i]   = dy_c_here_di*(1. - coef_up)   + (infos->H/2)*infos->dxdi[i]*(1.-h_j_up)  *df_plusdi;
            infos->dy_downdi[j][i] = dy_c_here_di*(1. - coef_down) + (infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }

        // DOWNNNNNNSTREAMMMM

        for (int i = i_trail; i < infos->M; i++){
            infos->y_up_t_plus_1[j][i]   =  (infos->H/2 + coef_coef*time_coef_t_plus_1)*coef_up   - coef_coef*time_coef_t_plus_1;
            infos->y_down_t_plus_1[j][i] =  (infos->H/2 - coef_coef*time_coef_t_plus_1)*coef_down - coef_coef*time_coef_t_plus_1;

            infos->dy_updj[j][i]   = (infos->H/2 + time_coef)*numerator_up   / denominator_up;
            infos->dy_downdj[j][i] = (infos->H/2 - time_coef)*numerator_down / denominator_down;
        }

        for (int i = 0; i < infos->M; i++){
            infos->dy_updt[j][i]   = (infos->y_up_t_plus_1[j][i]   - infos->y_up[j][i])   / infos->dt;
            infos->dy_downdt[j][i] = (infos->y_down_t_plus_1[j][i] - infos->y_down[j][i]) / infos->dt;

            infos->D_up[j][i]          = infos->dxdi[i]*infos->dy_updj[j][i];
            infos->D_down[j][i]        = infos->dxdi[i]*infos->dy_downdj[j][i];
            infos->D_up_t_plus_1[j][i]   = infos->D_up[j][i];
            infos->D_down_t_plus_1[j][i] = infos->D_down[j][i];
        }
    }
}

void update_y_grid_mmmmmh(struct grid_and_data *infos){

    int i_lead  = infos->i_lead;
    int i_trail = infos->i_trail;
    int di = 8;

    double time_coef          = fonction_time_1(infos, infos->t);
    double time_coef_t_plus_1 = fonction_time_1(infos, infos->t + infos->dt);

    for (int j = 0; j < infos->N; j++){

        double arg_up   = 1. - (double)j / (infos->N-1); 
        double arg_down = (double)j / (infos->N-1);       

        double coef_up   = 1. - tanh(infos->beta_plus *arg_up)   / tanh(infos->beta);
        double coef_down = tanh(infos->beta_minus*arg_down)  / tanh(infos->beta) -1.;

        double h_j_up   = coef_up;   
        double h_j_down = coef_down;

        double numerator_up   = infos->beta_plus  / (double)(infos->N-1);
        double denominator_up = tanh(infos->beta)*cosh(infos->beta_plus*arg_up)*cosh(infos->beta_plus*arg_up);
        double dh_updj_coef   = (infos->H/2)*numerator_up / denominator_up;

        double numerator_down   = infos->beta_minus / (double)(infos->N-1);
        double denominator_down = tanh(infos->beta)*cosh(infos->beta_minus*arg_down)*cosh(infos->beta_minus*arg_down);
        double dh_downdj_coef   = (infos->H/2)*numerator_down / denominator_down;


        for (int i = 0; i < infos->M; i++){
            infos->y_up[j][i]   = infos->y_up_t_plus_1[j][i];
            infos->y_down[j][i] = infos->y_down_t_plus_1[j][i];

            infos->dy_updi[j][i]   = 0.;
            infos->dy_downdi[j][i] = 0.;
        }

        // upstream 
        for (int i = 0; i < i_lead; i++){
            infos->y_up_t_plus_1[j][i]   = (infos->H/2 - time_coef_t_plus_1)*coef_up   + time_coef_t_plus_1;
            infos->y_down_t_plus_1[j][i] = (infos->H/2 + time_coef_t_plus_1)*coef_down + time_coef_t_plus_1;

            infos->dy_updj[j][i]   = (infos->H/2 - time_coef)*numerator_up   / denominator_up;
            infos->dy_downdj[j][i] = (infos->H/2 + time_coef)*numerator_down / denominator_down;
        }


        // Leading edge smoothing
        for (int i = i_lead; i < i_lead + di; i++){


            double fp = f_plus(infos, i);
            double fm = f_minus(infos, i);
            double df_plusdi  = -0.5*(infos->x[i]-infos->x_c)
                            / sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c));
            double df_minusdi = -df_plusdi;

            double step  = (double)(i - i_lead) / (double)di;
            double scale = step*step*(3. - 2.*step);  // 0 → 1

            // Valeur upstream (ce qu'on avait juste avant i_lead)
            double y_up_upstream   = (infos->H/2 - time_coef)*coef_up   + time_coef;
            double y_down_upstream = (infos->H/2 + time_coef)*coef_down + time_coef;

            // Valeur airfoil complet (ce qu'on aura à i_lead+di)
            double y_up_airfoil   = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *f_plus(infos, i);
            double y_down_airfoil = (infos->H/2)*coef_down + (infos->H/2)*arg_down*f_minus(infos, i);

            // Interpolation cubique entre les deux
            infos->y_up_t_plus_1[j][i]   = (1.-scale)*y_up_upstream   + scale*y_up_airfoil;
            infos->y_down_t_plus_1[j][i] = (1.-scale)*y_down_upstream + scale*y_down_airfoil;

            infos->dy_updj[j][i]   = dh_updj_coef   + scale*dh_updj_coef  *fp;
            infos->dy_downdj[j][i] = dh_downdj_coef - scale*dh_downdj_coef*fm;

            infos->dy_updi[j][i]   = scale*(infos->H/2)*infos->dxdi[i]*(1.-h_j_up)  *df_plusdi;
            infos->dy_downdi[j][i] = scale*(infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }

        // arifoil central
        for (int i = i_lead + di; i <= i_trail - di; i++){
            double df_plusdi  = -0.5*(infos->x[i]-infos->x_c)/ sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c));
            double df_minusdi = -df_plusdi;

            infos->y_up_t_plus_1[j][i]   = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *f_plus(infos, i);
            infos->y_down_t_plus_1[j][i] = (infos->H/2)*coef_down + (infos->H/2)*arg_down*f_minus(infos, i);;

            infos->dy_updj[j][i]   = dh_updj_coef   + dh_updj_coef  *f_plus(infos, i);
            infos->dy_downdj[j][i] = dh_downdj_coef - dh_downdj_coef*f_minus(infos, i);;

            infos->dy_updi[j][i]   = (infos->H/2)*infos->dxdi[i]*(1.-h_j_up)  *df_plusdi;
            infos->dy_downdi[j][i] = (infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }

        // Trailing edge smoothing  
        for (int i = i_trail - di; i <= i_trail; i++){


            double fp = f_plus(infos, i);
            double fm = f_minus(infos, i);
            double df_plusdi  = -0.5*(infos->x[i]-infos->x_c)
                            / sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c));
            double df_minusdi = -df_plusdi;

            double step  = (double)(i_trail - i) / (double)di;
            double scale = step*step*(3. - 2.*step);  // 1 → 0



            // Valeur downstream (ce qu'on aura juste après i_trail)
            double y_up_downstream   = (infos->H/2 + time_coef)*coef_up   - time_coef;
            double y_down_downstream = (infos->H/2 - time_coef)*coef_down - time_coef;

            // Valeur airfoil complet (ce qu'on avait à i_trail-di)
            double y_up_airfoil   = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *f_plus(infos, i);
            double y_down_airfoil = (infos->H/2)*coef_down + (infos->H/2)*arg_down*f_minus(infos, i);

            // Interpolation cubique entre airfoil et downstream
            infos->y_up_t_plus_1[j][i]   = scale*y_up_airfoil   + (1.-scale)*y_up_downstream;
            infos->y_down_t_plus_1[j][i] = scale*y_down_airfoil + (1.-scale)*y_down_downstream;



            infos->dy_updj[j][i]   = dh_updj_coef   + scale*dh_updj_coef  *fp;
            infos->dy_downdj[j][i] = dh_downdj_coef - scale*dh_downdj_coef*fm;

            infos->dy_updi[j][i]   = scale*(infos->H/2)*infos->dxdi[i]*(1.-h_j_up)  *df_plusdi;
            infos->dy_downdi[j][i] = scale*(infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }
        // downstream
        for (int i = i_trail+1; i < infos->M; i++){
            infos->y_up_t_plus_1[j][i]   =  (infos->H/2 + time_coef_t_plus_1)*coef_up   - time_coef_t_plus_1;
            infos->y_down_t_plus_1[j][i] =  (infos->H/2 - time_coef_t_plus_1)*coef_down - time_coef_t_plus_1;

            infos->dy_updj[j][i]   = (infos->H/2 + time_coef)*numerator_up   / denominator_up;
            infos->dy_downdj[j][i] = (infos->H/2 - time_coef)*numerator_down / denominator_down;
        }

        for (int i = 0; i < infos->M; i++){
            infos->dy_updt[j][i]   = (infos->y_up_t_plus_1[j][i]   - infos->y_up[j][i])   / infos->dt;
            infos->dy_downdt[j][i] = (infos->y_down_t_plus_1[j][i] - infos->y_down[j][i]) / infos->dt;

            infos->D_up[j][i]          = infos->dxdi[i]*infos->dy_updj[j][i];
            infos->D_down[j][i]        = infos->dxdi[i]*infos->dy_downdj[j][i];
            infos->D_up_t_plus_1[j][i]   = infos->D_up[j][i];
            infos->D_down_t_plus_1[j][i] = infos->D_down[j][i];
        }
    }
}

void update_y_grid(struct grid_and_data *infos){

    int i_lead  = infos->i_lead;
    int i_trail = infos->i_trail;
    int di = 8;

    double time_coef          = fonction_time_1(infos, infos->t);
    double time_coef_t_plus_1 = fonction_time_1(infos, infos->t + infos->dt);

    for (int j = 0; j < infos->N; j++){

        double arg_up   = 1. - (double)j / (infos->N-1);
        double arg_down = (double)j / (infos->N-1);

        double coef_up   = 1. - tanh(infos->beta_plus *arg_up)   / tanh(infos->beta);
        double coef_down = tanh(infos->beta_minus*arg_down) / tanh(infos->beta) - 1.;

        double h_j_up   = coef_up;
        double h_j_down = coef_down;

        double numerator_up   = infos->beta_plus  / (double)(infos->N-1);
        double denominator_up = tanh(infos->beta)*cosh(infos->beta_plus*arg_up)*cosh(infos->beta_plus*arg_up);
        double dh_updj_coef   = (infos->H/2)*numerator_up / denominator_up;

        double numerator_down   = infos->beta_minus / (double)(infos->N-1);
        double denominator_down = tanh(infos->beta)*cosh(infos->beta_minus*arg_down)*cosh(infos->beta_minus*arg_down);
        double dh_downdj_coef   = (infos->H/2)*numerator_down / denominator_down;

        for (int i = 0; i < infos->M; i++){
            infos->y_up[j][i]   = infos->y_up_t_plus_1[j][i];
            infos->y_down[j][i] = infos->y_down_t_plus_1[j][i];
            infos->dy_updi[j][i]   = 0.;
            infos->dy_downdi[j][i] = 0.;
        }

        // upstream
        for (int i = 0; i < i_lead; i++){
            infos->y_up_t_plus_1[j][i]   = (infos->H/2 - time_coef_t_plus_1)*coef_up   + time_coef_t_plus_1;
            infos->y_down_t_plus_1[j][i] = (infos->H/2 + time_coef_t_plus_1)*coef_down + time_coef_t_plus_1;

            infos->dy_updj[j][i]   = (infos->H/2 - time_coef)*numerator_up   / denominator_up;
            infos->dy_downdj[j][i] = (infos->H/2 + time_coef)*numerator_down / denominator_down;
        }

        // Leading edge smoothing
        for (int i = i_lead; i < i_lead + di; i++){
            double step  = (double)(i - i_lead) / (double)di;
            double scale = step*step*(3. - 2.*step);  // 0 → 1

            double fp = f_plus(infos, i);
            double fm = f_minus(infos, i);
            double df_plusdi  = -0.5*(infos->x[i]-infos->x_c)
                               / sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c));
            double df_minusdi = -df_plusdi;

            // y_t_plus_1 : interpolation avec time_coef_t_plus_1
            double y_up_upstream_tp1   = (infos->H/2 - time_coef_t_plus_1)*coef_up   + time_coef_t_plus_1;
            double y_down_upstream_tp1 = (infos->H/2 + time_coef_t_plus_1)*coef_down + time_coef_t_plus_1;

            double y_up_airfoil        = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *fp;
            double y_down_airfoil      = (infos->H/2)*coef_down + (infos->H/2)*arg_down*fm;

            infos->y_up_t_plus_1[j][i]   = (1.-scale)*y_up_upstream_tp1   + scale*y_up_airfoil;
            infos->y_down_t_plus_1[j][i] = (1.-scale)*y_down_upstream_tp1 + scale*y_down_airfoil;

            // dérivées : interpolation avec time_coef
            infos->dy_updj[j][i]   = (1.-scale)*(infos->H/2 - time_coef)*numerator_up  /denominator_up
                                    + scale*(dh_updj_coef   + dh_updj_coef  *fp);
            infos->dy_downdj[j][i] = (1.-scale)*(infos->H/2 + time_coef)*numerator_down/denominator_down
                                    + scale*(dh_downdj_coef - dh_downdj_coef*fm);

            infos->dy_updi[j][i]   = scale*(infos->H/2)*infos->dxdi[i]*(1.-h_j_up)  *df_plusdi;
            infos->dy_downdi[j][i] = scale*(infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }

        // airfoil central
        for (int i = i_lead + di; i <= i_trail - di; i++){
            double df_plusdi  = -0.5*(infos->x[i]-infos->x_c)
                               / sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c));
            double df_minusdi = -df_plusdi;

            infos->y_up_t_plus_1[j][i]   = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *f_plus(infos, i);
            infos->y_down_t_plus_1[j][i] = (infos->H/2)*coef_down + (infos->H/2)*arg_down*f_minus(infos, i);

            infos->dy_updj[j][i]   = dh_updj_coef   + dh_updj_coef  *f_plus(infos, i);
            infos->dy_downdj[j][i] = dh_downdj_coef - dh_downdj_coef*f_minus(infos, i);

            infos->dy_updi[j][i]   = (infos->H/2)*infos->dxdi[i]*(1.-h_j_up)  *df_plusdi;
            infos->dy_downdi[j][i] = (infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }

        // Trailing edge smoothing
        for (int i = i_trail - di; i <= i_trail; i++){
            double step  = (double)(i_trail - i) / (double)di;
            double scale = step*step*(3. - 2.*step);  // 1 → 0

            double fp = f_plus(infos, i);
            double fm = f_minus(infos, i);
            double df_plusdi  = -0.5*(infos->x[i]-infos->x_c)
                               / sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c));
            double df_minusdi = -df_plusdi;

            // y_t_plus_1 : interpolation avec time_coef_t_plus_1
            double y_up_downstream_tp1   = (infos->H/2 + time_coef_t_plus_1)*coef_up   - time_coef_t_plus_1;
            double y_down_downstream_tp1 = (infos->H/2 - time_coef_t_plus_1)*coef_down - time_coef_t_plus_1;
            double y_up_airfoil          = (infos->H/2)*coef_up   + (infos->H/2)*arg_up  *fp;
            double y_down_airfoil        = (infos->H/2)*coef_down + (infos->H/2)*arg_down*fm;

            infos->y_up_t_plus_1[j][i]   = scale*y_up_airfoil   + (1.-scale)*y_up_downstream_tp1;
            infos->y_down_t_plus_1[j][i] = scale*y_down_airfoil + (1.-scale)*y_down_downstream_tp1;

            // dérivées : interpolation avec time_coef
            infos->dy_updj[j][i]   = scale*(dh_updj_coef   + dh_updj_coef  *fp)
                                    + (1.-scale)*(infos->H/2 + time_coef)*numerator_up  /denominator_up;
            infos->dy_downdj[j][i] = scale*(dh_downdj_coef - dh_downdj_coef*fm)
                                    + (1.-scale)*(infos->H/2 - time_coef)*numerator_down/denominator_down;

            infos->dy_updi[j][i]   = scale*(infos->H/2)*infos->dxdi[i]*(1.-h_j_up)  *df_plusdi;
            infos->dy_downdi[j][i] = scale*(infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }

        // downstream
        for (int i = i_trail+1; i < infos->M; i++){
            infos->y_up_t_plus_1[j][i]   = (infos->H/2 + time_coef_t_plus_1)*coef_up   - time_coef_t_plus_1;
            infos->y_down_t_plus_1[j][i] = (infos->H/2 - time_coef_t_plus_1)*coef_down - time_coef_t_plus_1;

            infos->dy_updj[j][i]   = (infos->H/2 + time_coef)*numerator_up   / denominator_up;
            infos->dy_downdj[j][i] = (infos->H/2 - time_coef)*numerator_down / denominator_down;
        }

        // dérivées temporelles et Jacobiens
        for (int i = 0; i < infos->M; i++){
            infos->dy_updt[j][i]   = (infos->y_up_t_plus_1[j][i]   - infos->y_up[j][i])   / infos->dt;
            infos->dy_downdt[j][i] = (infos->y_down_t_plus_1[j][i] - infos->y_down[j][i]) / infos->dt;

            infos->D_up[j][i]            = infos->dxdi[i]*infos->dy_updj[j][i];
            infos->D_down[j][i]          = infos->dxdi[i]*infos->dy_downdj[j][i];
            infos->D_up_t_plus_1[j][i]   = infos->D_up[j][i];
            infos->D_down_t_plus_1[j][i] = infos->D_down[j][i];
        }
    }
}

struct grid_and_data init_grid_and_data(double M_inf, double grow_rate_M_inf, double t_end ,int step_max, double amplitude_oscillation, double frequence_oscillation, double angle_of_attack){

    /*
        pre:
         M_inf : the Mach number of the flow must be positive,
         t_end : the simulation must run for a positive amount of time,
         step_max > 0 => after how many steps we write the results on the files, usually after t_end*1000, steps, it produces an animation of 3 secondes
         amplitude_oscillation:  the amplitude of the oscillation must be non-negative, the amplitudes
         frequence_oscillation >= 0, 
         angle_of_attack >= 0
        post: returns a struct grid_and_data with the physical parameters, geometric parameters, grid and derivatives initialized based on the input parameters and some predefined values.
    */

    struct grid_and_data infos_local;
    struct grid_and_data *infos = &infos_local;

    infos->M_inf_init = M_inf;
    infos->M_inf = M_inf;
    infos->grow_rate_M_inf = grow_rate_M_inf;

    // paramètres temporels
    infos->t = 0.;
    infos->dt = 1e-5; 
    infos->t_end = t_end;
    infos->step_max = step_max;


    infos->C = .1;
    infos->eps_min = 2.; 

    infos->gamma = 1.4;
    infos->T_inf = 300.0;
    infos->p_inf = 1e5;
    infos->R_fluid = 287.1;

    infos->rho_inf =infos->p_inf/ (infos->R_fluid *infos->T_inf);
    infos->c_sound_inf = sqrt(infos->gamma*infos->R_fluid*infos->T_inf);

    infos->T_0_inf = infos->T_inf*(1.+0.5*(infos->gamma-1.)*infos->M_inf*infos->M_inf);
    infos->p_0_inf = infos->p_inf*pow(infos->T_0_inf / infos->T_inf, infos->gamma/(infos->gamma-1) );

    // paramètres spaciaux de la grille

    double c = 1.;
    double L = 5.;

    infos->c = c;
    infos->R = 3*c;
    infos->H = 4*c;
    infos->y_c = 0.0;
    infos->x_c = c/2;

    infos->L = L;
    infos->L_u = 1.5;
    


    infos->M = 150; 
    infos->i_1 = (int)infos->M/5;
    infos->i_2 = (int)infos->M*4/5;
    infos->i_lead = 0;// déterminé  dans une autre fct
    infos->i_trail = 0; // déterminé  dans une autre fct
    infos->N = 30;  

    infos->beta_1 = 0.81106561; // déterminé va statement & python
    infos->beta_3 = 1.44500835; // déterminé va statement & python 

    /*

    infos->M = 450; 
    infos->i_1 = (int)infos->M/5;
    infos->i_2 = (int)infos->M*4/5;
    infos->i_lead = 0;// déterminé  dans une autre fct
    infos->i_trail = 0; // déterminé  dans une autre fct
    infos->N = 90;  

    infos->beta_1 = 0.81106561; // déterminé va statement & python
    infos->beta_3 = 1.42774587; // déterminé va statement & python  
    infos->dt = 5e-6;
    // laissés tels quel ça change rien
    */

    infos->beta_plus = 1.; 
    infos->beta_minus = 1.;
    infos->beta = 1.;

    // paramètres de l'aile
    infos->AOA = angle_of_attack*M_PI/180.; // passe de deg à rad
    infos->mu = amplitude_oscillation/2.;
    infos->omega = frequence_oscillation; 

    // grilles ,dérivées et caractéristiques
    infos->x = (double *)calloc(infos->M, sizeof(double));
    infos->y_up = (double **)calloc(infos->N, sizeof(double*));
    infos->y_down = (double **)calloc(infos->N, sizeof(double*));

    infos->y_up_t_plus_1 = (double **)calloc(infos->N, sizeof(double*));
    infos->y_down_t_plus_1 = (double **)calloc(infos->N, sizeof(double*));
    infos->dxdi = (double *)calloc(infos->M, sizeof(double));

    infos->dy_updi = (double **)calloc(infos->N, sizeof(double*));
    infos->dy_updj = (double **)calloc(infos->N, sizeof(double*));
    infos->dy_updt = (double **)calloc(infos->N, sizeof(double*));

    infos->dy_downdi = (double **)calloc(infos->N, sizeof(double*));
    infos->dy_downdj = (double **)calloc(infos->N, sizeof(double*));
    infos->dy_downdt = (double **)calloc(infos->N, sizeof(double*));

    infos->D_up = (double **)calloc(infos->N, sizeof(double*));
    infos->D_down = (double **)calloc(infos->N, sizeof(double*));
    infos->D_up_t_plus_1 = (double **)calloc(infos->N, sizeof(double*));
    infos->D_down_t_plus_1 = (double **)calloc(infos->N, sizeof(double*));

    infos->s_up = (double***)calloc(4, sizeof(double**));
    infos->s_down = (double***)calloc(4, sizeof(double**));

    infos->s_star_up = (double***)calloc(4, sizeof(double**));
    infos->s_star_down = (double***)calloc(4, sizeof(double**));

    infos->F_x_up = (double***)calloc(4, sizeof(double**));
    infos->F_x_down = (double***)calloc(4, sizeof(double**));
    infos->F_x_star_up = (double***)calloc(4, sizeof(double**));
    infos->F_x_star_down = (double***)calloc(4, sizeof(double**));
    
    infos->F_y_up = (double***)calloc(4, sizeof(double**));
    infos->F_y_down = (double***)calloc(4, sizeof(double**));
    infos->F_y_star_up = (double***)calloc(4, sizeof(double**));
    infos->F_y_star_down = (double***)calloc(4, sizeof(double**));

    infos->F_i_up = (double***)calloc(4, sizeof(double**));
    infos->F_i_down = (double***)calloc(4, sizeof(double**));
    infos->F_i_star_up = (double***)calloc(4, sizeof(double**));
    infos->F_i_star_down = (double***)calloc(4, sizeof(double**));

    infos->F_j_up = (double***)calloc(4, sizeof(double**));
    infos->F_j_down = (double***)calloc(4, sizeof(double**));
    infos->F_j_star_up = (double***)calloc(4, sizeof(double**));
    infos->F_j_star_down = (double***)calloc(4, sizeof(double**));

    for (int i = 0; i < 4; i++) {
        infos->s_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->s_down[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->s_star_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->s_star_down[i] = (double **)calloc(infos->N, sizeof(double*));

        infos->F_x_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_x_down[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_x_star_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_x_star_down[i] = (double **)calloc(infos->N, sizeof(double*));

        infos->F_y_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_y_down[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_y_star_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_y_star_down[i] = (double **)calloc(infos->N, sizeof(double*));

        infos->F_i_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_i_down[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_i_star_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_i_star_down[i] = (double **)calloc(infos->N, sizeof(double*));

        infos->F_j_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_j_down[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_j_star_up[i] = (double **)calloc(infos->N, sizeof(double*));
        infos->F_j_star_down[i] = (double **)calloc(infos->N, sizeof(double*));
    }

    for (int j =0; j <infos->N;j++){
        infos->y_up[j] = (double*)calloc(infos->M, sizeof(double));
        infos->y_down[j] = (double*)calloc(infos->M, sizeof(double));

        infos->y_up_t_plus_1[j] = (double*)calloc(infos->M, sizeof(double));
        infos->y_down_t_plus_1[j] = (double*)calloc(infos->M, sizeof(double));

        infos->D_up[j] = (double*)calloc(infos->M, sizeof(double));
        infos->D_down[j] = (double*)calloc(infos->M, sizeof(double));
        infos->D_up_t_plus_1[j] = (double*)calloc(infos->M, sizeof(double));
        infos->D_down_t_plus_1[j] = (double*)calloc(infos->M, sizeof(double));

        infos->dy_updi[j] = (double*)calloc(infos->M, sizeof(double));
        infos->dy_updj[j] = (double*)calloc(infos->M, sizeof(double));
        infos->dy_updt[j] = (double*)calloc(infos->M, sizeof(double));

        infos->dy_downdi[j] = (double*)calloc(infos->M, sizeof(double));
        infos->dy_downdj[j] = (double*)calloc(infos->M, sizeof(double));
        infos->dy_downdt[j] = (double*)calloc(infos->M, sizeof(double));

        for (int i =0; i<4;i++){
            infos->s_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->s_down[i][j] = (double*)calloc(infos->M, sizeof(double));

            infos->s_star_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->s_star_down[i][j] = (double*)calloc(infos->M, sizeof(double));

            infos->F_x_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_x_down[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_x_star_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_x_star_down[i][j] = (double*)calloc(infos->M, sizeof(double));

            infos->F_y_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_y_down[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_y_star_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_y_star_down[i][j] = (double*)calloc(infos->M, sizeof(double));

            infos->F_i_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_i_down[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_i_star_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_i_star_down[i][j] = (double*)calloc(infos->M, sizeof(double));

            infos->F_j_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_j_down[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_j_star_up[i][j] = (double*)calloc(infos->M, sizeof(double));
            infos->F_j_star_down[i][j] = (double*)calloc(infos->M, sizeof(double));
        }
    }

    
    return *infos;
}

void free_data(struct grid_and_data *infos){
    free(infos->x);
    free(infos->dxdi);

    for (int j = 0; j < infos->N; j++) {
        free(infos->y_up[j]);
        free(infos->y_down[j]);

        free(infos->y_up_t_plus_1[j]);
        free(infos->y_down_t_plus_1[j]);

        free(infos->D_up[j]);
        free(infos->D_down[j]);
        free(infos->D_up_t_plus_1[j]);
        free(infos->D_down_t_plus_1[j]);

        free(infos->dy_updi[j]);
        free(infos->dy_updj[j]);
        free(infos->dy_updt[j]);

        free(infos->dy_downdi[j]);
        free(infos->dy_downdj[j]);
        free(infos->dy_downdt[j]);



        for (int i = 0; i < 4; i++){
            free(infos->s_up[i][j]);
            free(infos->s_down[i][j]);
            free(infos->s_star_up[i][j]);
            free(infos->s_star_down[i][j]);

            free(infos->F_x_up[i][j]);
            free(infos->F_x_down[i][j]);
            free(infos->F_x_star_up[i][j]);
            free(infos->F_x_star_down[i][j]);

            free(infos->F_y_up[i][j]);
            free(infos->F_y_down[i][j]);
            free(infos->F_y_star_up[i][j]);
            free(infos->F_y_star_down[i][j]);

            free(infos->F_i_up[i][j]);
            free(infos->F_i_down[i][j]);
            free(infos->F_i_star_up[i][j]);
            free(infos->F_i_star_down[i][j]);

            free(infos->F_j_up[i][j]);
            free(infos->F_j_down[i][j]);
            free(infos->F_j_star_up[i][j]);
            free(infos->F_j_star_down[i][j]);
        }
    }

    for (int i = 0; i < 4; i++) {
        free(infos->s_up[i]);
        free(infos->s_down[i]);
        free(infos->s_star_up[i]);
        free(infos->s_star_down[i]);

        free(infos->F_x_up[i]);
        free(infos->F_x_down[i]);
        free(infos->F_x_star_up[i]);
        free(infos->F_x_star_down[i]);

        free(infos->F_y_up[i]);
        free(infos->F_y_down[i]);
        free(infos->F_y_star_up[i]);
        free(infos->F_y_star_down[i]);

        free(infos->F_i_up[i]);
        free(infos->F_i_down[i]);
        free(infos->F_i_star_up[i]);
        free(infos->F_i_star_down[i]);

        free(infos->F_j_up[i]);
        free(infos->F_j_down[i]);
        free(infos->F_j_star_up[i]);
        free(infos->F_j_star_down[i]);
    }

    free(infos->y_up);
    free(infos->y_down);
    free(infos->y_up_t_plus_1);
    free(infos->y_down_t_plus_1);

    free(infos->D_up);
    free(infos->D_down);
    free(infos->D_up_t_plus_1);
    free(infos->D_down_t_plus_1);


    free(infos->dy_updi);
    free(infos->dy_updj);
    free(infos->dy_updt);

    free(infos->dy_downdi);
    free(infos->dy_downdj);
    free(infos->dy_downdt);

    free(infos->s_up);
    free(infos->s_down);
    free(infos->s_star_up);
    free(infos->s_star_down);

    free(infos->F_x_up);
    free(infos->F_x_down);
    free(infos->F_x_star_up);
    free(infos->F_x_star_down);

    free(infos->F_y_up);
    free(infos->F_y_down);
    free(infos->F_y_star_up);
    free(infos->F_y_star_down);

    free(infos->F_i_up);
    free(infos->F_i_down);
    free(infos->F_i_star_up);
    free(infos->F_i_star_down);

    free(infos->F_j_up);
    free(infos->F_j_down);
    free(infos->F_j_star_up);
    free(infos->F_j_star_down);
}

#endif 