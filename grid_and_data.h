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
    bool hdd_storage;
    double resolution;
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
    double number_of_wrotten_steps; // 1000 steps pr une anumation de 33 secondes pour 1000 steps, 100 steps pr une animation de 3 secondes
    double write_interval_steps;

    // some parameter we must choose carefuly
    double beta_1;
    double beta_3;
    double beta_minus;
    double beta_plus;
    double beta;

    double AOA;
    double AOA_current;
    double t_ramp;
    double AOA_grow_rate;
    double mu;
    double omega;

    // grids
    double* x;  
    double* y_up;
    double* y_down;


    //derivatives and other stuff

    double* D_up;
    double* D_down;
    double* D_up_t_plus_1; // à changer quand ça bougera au cours du temps
    double* D_down_t_plus_1; // à changer quand ça bougera au cours du temps

    double* dxdi;
    double* dydi_up;
    double* dydj_up;
    double* dydt_up;
    double* dydi_down;
    double* dydj_down;
    double* dydt_down;
    //no dxdj and no dxdt


    // physical dimensions in terms of s and F
    
    double *s_up;
    double *s_down;
    double *s_star_up;
    double *s_star_down;

    double *F_i_up;
    double *F_i_down;
    double *F_i_star_up; 
    double *F_i_star_down; 

    double *F_j_up;
    double *F_j_down;
    double *F_j_star_up; 
    double *F_j_star_down; 
};
  

void leading_and_trailing_edge_finder(struct grid_and_data *infos){
    /*
        pre :
            infos : the structure we use
        post :
            returns the x  indexes of the leading and trailing edges
    */
    int i = 0;

    while (infos->x[i] < 0.){
        
        i +=1;
    }
    infos->i_lead = i;
    while (infos->x[i] <1.){
        
        i +=1;
    }
    infos->i_trail = i;

}

// case 1

double f_plus_case1(struct grid_and_data*infos, double x_i, double t){
    /*pre :
            infos : the structure we use
            x_i : the x_index where we compute the f_plus value
            t : the time at which  we compute the f_plus value
        post :
            returns the y value of the point x_i at time t. it represents our upper part of the airfoil
    */
    double value_space = sqrt(-(x_i-infos->x_c)*(x_i-infos->x_c)+infos->R*infos->R)- sqrt(-infos->x_c*infos->x_c+infos->R*infos->R) ;
    double value_AoA = (1.-2.*x_i/infos->c)*infos->c*tan(infos->AOA_current);
    double value_time = (1.-2.*x_i/infos->c)* infos->mu*infos->c*sin(infos->omega*t) ;
    return (value_space+value_AoA+ value_time)  ;
}

double f_minus_case1(struct grid_and_data*infos, double x_i, double t){
    /*
        pre :
            infos : the structure we use
            x_i : the x_index where we compute the f_minus value
            t : the time at which  we compute the f_minus value
        post :
            returns the y value of the point x_i at time t. it represents our lower part of the airfoil
    */
    double value_space = sqrt(-(x_i-infos->x_c)*(x_i-infos->x_c)+infos->R*infos->R)- sqrt(-infos->x_c*infos->x_c+infos->R*infos->R);
    double value_AoA = (1.-2.*x_i/infos->c)*infos->c*tan(infos->AOA_current);
    double value_time = (1.-2.*x_i/infos->c)* infos->mu*infos->c*sin(infos->omega*t) ;
    return ( -value_space+value_AoA+value_time);
}

double ddt_f_plus_case1(struct grid_and_data*infos, double x_i){
    /*
        pre :
            infos : the structure we use
            x_i : the x_index where we compute the f_plus value
           
        post :
            returns the y value of the point x_i a. it represents the time derivative of the upper part of the airfoil
    */
    double value_space = 0.;
    double value_AoA =  (1.-2.*x_i/infos->c) * infos->c *infos->AOA_grow_rate/(cos(infos->AOA_current)*cos(infos->AOA_current));
    double value_time = (1.-2.*x_i/infos->c)*infos->omega* infos->mu*infos->c*cos(infos->omega*infos->t) ;
    return (value_space+value_AoA+value_time)  ;
}

double ddt_f_minus_case1(struct grid_and_data*infos, double x_i){
    /*
        pre :
            infos : the structure we use
            x_i : the x_index where we compute the f_plus value
           
        post :
            returns the y value of the point x_i a. it represents the time derivative of the lower part of the airfoil
    */
    double value_space = 0.;
    double value_AoA =  (1.-2.*x_i/infos->c) * infos->c *infos->AOA_grow_rate/(cos(infos->AOA_current)*cos(infos->AOA_current));
    double value_time = (1.-2.*x_i/infos->c)*infos->omega* infos->mu*infos->c*cos(infos->omega*infos->t) ;
    return (value_space+value_AoA+value_time)  ;
}

double ddx_f_plus_case1(struct grid_and_data*infos, double x_i){
    /*
        pre :
            infos : the structure we use
            x_i : the x_index where we compute the f_plus value
           
        post :
            returns the y value of the point x_i a. it represent the space derivative of the upper part of the airfoil
    */
    double value_space = -(x_i-infos->x_c )/sqrt(-(x_i-infos->x_c)*(x_i-infos->x_c)+infos->R*infos->R);
    double value_AoA = -2.*tan(infos->AOA_current);
    double value_time = (-2./infos->c)* infos->mu*infos->c*sin(infos->omega*infos->t) ;
    return (value_space+value_AoA+value_time)  ;
}

double ddx_f_minus_case1(struct grid_and_data*infos, double x_i){
     /*
        pre :
            infos : the structure we use
            x_i : the x_index where we compute the f_plus value
        post :
            returns the y value of the point x_i a. it represents the space derivative of the lower part of the airfoil
    */
    double value_space = -(x_i-infos->x_c )/sqrt(-(x_i-infos->x_c)*(x_i-infos->x_c)+infos->R*infos->R);
    double value_AoA = -2.*tan(infos->AOA_current);
    double value_time = (-2./infos->c)* infos->mu*infos->c*sin(infos->omega*infos->t) ;

    return (-value_space+value_AoA+value_time)  ;
}

// case 2 stretching arfoil
///TDÖKJIUVUC

// Grid stuff

void set_x_grid(struct grid_and_data *infos){ 
     /*
        pre :
            infos : the structure we use
        post :
            put the x and dxdi values of the grid in infos. x is only one dimesion
    */
    
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

void set_y_grid(struct grid_and_data *infos){ 
    /* xpre :
            infos : the structure we use
        post :
            put the y values of the grid in the upper and lower part. it also compute the derivatives of yin infos.
    */

    //on calcule les points en y, on calcule les dérivées dans update_y_grid
    int M = infos->M;
    int N = infos->N;
    int i_lead  = infos->i_lead;
    int i_trail = infos->i_trail;

    for (int j = 0; j < infos->N; j++){

        double H = infos->H;

        double arg_up   = 1. - (double)j / (double)(infos->N-1.); 
        double arg_down = (double)j / (double)(infos->N-1.);       

        double h_up   = 1. - tanh(infos->beta_plus *arg_up)   / tanh(infos->beta);
        double h_down = tanh(infos->beta_minus*arg_down)  / tanh(infos->beta) -1.;
        
        double dhdj_den  =  ((double)infos->N-1.)*tanh(infos->beta);

        double dhdj_up   =  infos->beta_plus / (dhdj_den*cosh(infos->beta_plus *arg_up)*cosh(infos->beta_plus *arg_up));
        double dhdj_down =  infos->beta_minus/ (dhdj_den*cosh(infos->beta_minus *arg_down)*cosh(infos->beta_minus *arg_down));
        
        for (int i = 0; i < i_lead; i++){
            double x_i = 0.;
            infos->y_up[j*M+i] = 0.5*H*h_up + arg_up*f_plus_case1(infos, x_i, infos->t);
            infos->y_down[j*M+i] = 0.5*H*h_down + arg_down*f_minus_case1(infos, x_i, infos->t);

            infos->dydt_up[j*M+i] = (1. - h_up)*ddt_f_plus_case1(infos, x_i);
            infos->dydt_down[j*M+i] = (1. + h_down)*ddt_f_minus_case1(infos, x_i);

            infos->dydi_up[j*M+i]   = 0.;
            infos->dydi_down[j*M+i] = 0.;

            infos->dydj_up[j*M+i] = dhdj_up*(0.5*H - f_plus_case1(infos, x_i, infos->t));
            infos->dydj_down[j*M+i] = dhdj_down*(0.5*H + f_minus_case1(infos, x_i, infos->t));

            infos->D_up[j*M+i]   = infos->dxdi[i]*infos->dydj_up[j*M+i];
            infos->D_down[j*M+i] = infos->dxdi[i]*infos->dydj_down[j*M+i];

            infos->D_up_t_plus_1[j*M+i] = infos->dxdi[i]*dhdj_up*(0.5*H - f_plus_case1(infos, x_i, infos->t+infos->dt));
            infos->D_down_t_plus_1[j*M+i] = infos->dxdi[i]*dhdj_down*(0.5*H + f_minus_case1(infos, x_i, infos->t+infos->dt));
        }
        

        // arifoil central
        for (int i = i_lead; i < i_trail; i++){
            double x_i = infos->x[i];
            infos->y_up[j*M+i] = 0.5*H*h_up + arg_up*f_plus_case1(infos, x_i, infos->t);
            infos->y_down[j*M+i] = 0.5*H*h_down + arg_down*f_minus_case1(infos, x_i, infos->t);

            infos->dydt_up[j*M+i] = (1. - h_up)*ddt_f_plus_case1(infos, x_i);
            infos->dydt_down[j*M+i] = (1. + h_down)*ddt_f_minus_case1(infos, x_i);

            infos->dydi_up[j*M+i]   = infos->dxdi[i]*( 1. - h_up)*ddx_f_plus_case1(infos, x_i);
            infos->dydi_down[j*M+i] = infos->dxdi[i]*( 1. + h_down)*ddx_f_minus_case1(infos, x_i);

            infos->dydj_up[j*M+i] = dhdj_up*(0.5*H - f_plus_case1(infos, x_i, infos->t));
            infos->dydj_down[j*M+i] = dhdj_down*(0.5*H + f_minus_case1(infos, x_i, infos->t));

            infos->D_up[j*M+i]   = infos->dxdi[i]*infos->dydj_up[j*M+i];
            infos->D_down[j*M+i] = infos->dxdi[i]*infos->dydj_down[j*M+i];

            infos->D_up_t_plus_1[j*M+i] = infos->dxdi[i]*dhdj_up*(0.5*H - f_plus_case1(infos, x_i, infos->t+infos->dt));
            infos->D_down_t_plus_1[j*M+i] = infos->dxdi[i]*dhdj_down*(0.5*H + f_minus_case1(infos, x_i, infos->t+infos->dt));
        }

        // downstream
        for (int i = i_trail; i < infos->M; i++){
            double x_i = 1.;
            infos->y_up[j*M+i] = 0.5*H*h_up + arg_up*f_plus_case1(infos, x_i, infos->t);
            infos->y_down[j*M+i] = 0.5*H*h_down + arg_down*f_minus_case1(infos, x_i, infos->t);

            infos->dydt_up[j*M+i] = (1. - h_up)*ddt_f_plus_case1(infos, x_i);
            infos->dydt_down[j*M+i] = (1. + h_down)*ddt_f_minus_case1(infos, x_i);

            infos->dydi_up[j*M+i]   = 0.;
            infos->dydi_down[j*M+i] = 0.;

            infos->dydj_up[j*M+i] = dhdj_up*(0.5*H - f_plus_case1(infos, x_i, infos->t));
            infos->dydj_down[j*M+i] = dhdj_down*(0.5*H + f_minus_case1(infos, x_i, infos->t));

            infos->D_up[j*M+i]   = infos->dxdi[i]*infos->dydj_up[j*M+i];
            infos->D_down[j*M+i] = infos->dxdi[i]*infos->dydj_down[j*M+i];

            infos->D_up_t_plus_1[j*M+i] = infos->dxdi[i]*dhdj_up*(0.5*H - f_plus_case1(infos, x_i, infos->t+infos->dt));
            infos->D_down_t_plus_1[j*M+i] = infos->dxdi[i]*dhdj_down*(0.5*H + f_minus_case1(infos, x_i, infos->t+infos->dt));
        }
    }


    
    int j_up = 0;
    int j_down = infos->N-1;

    int di_lead = 4;
    int di = 4;

    double H = infos->H;

    double arg_up   = 1. - (double)j_up / (double)(infos->N-1.); 
    double arg_down = (double)j_down / (double)(infos->N-1.);       

    double h_up   = 1. - tanh(infos->beta_plus *arg_up)   / tanh(infos->beta);
    double h_down = tanh(infos->beta_minus*arg_down)  / tanh(infos->beta) -1.;
    
    double dhdj_den  =  ((double)infos->N-1.)*tanh(infos->beta);

    double dhdj_up   =  infos->beta_plus / (dhdj_den*cosh(infos->beta_plus *arg_up)*cosh(infos->beta_plus *arg_up));
    double dhdj_down =  infos->beta_minus/ (dhdj_den*cosh(infos->beta_minus *arg_down)*cosh(infos->beta_minus *arg_down));
    
   
    // smoothing at the two bounda ry rows only (leading and trailing)
    // leading smoothing (only j_up row)
    for (int i = i_lead; i < i_lead + di_lead; ++i){
        double s = (double)(i - i_lead) / (double)di_lead;
        double alpha = s*s*(3.0 - 2.0*s); // smoothstep
        double x_i = infos->x[i];
        double w0 = 1.0 - alpha; // weight for endpoint value (x=0)
        double wx = alpha;       // weight for airfoil value (x_i)

        infos->y_up[j_up*M + i]   = 0.5*H*h_up + arg_up*( w0 * f_plus_case1(infos, 0.0, infos->t) + wx * f_plus_case1(infos, x_i, infos->t) );
        infos->y_down[j_down*M+i] = 0.5*H*h_down + arg_down*( w0 * f_minus_case1(infos, 0.0, infos->t) + wx * f_minus_case1(infos, x_i, infos->t) );

        infos->dydt_up[j_up*M + i]   = (1.0 - h_up)*( w0 * ddt_f_plus_case1(infos, 0.0) + wx * ddt_f_plus_case1(infos, x_i) );
        infos->dydt_down[j_down*M + i] = (1.0 + h_down)*( w0 * ddt_f_minus_case1(infos, 0.0) + wx * ddt_f_minus_case1(infos, x_i) );

        infos->dydi_up[j_up*M + i]   = infos->dxdi[i] * (1.0 - h_up) * ( w0 * ddx_f_plus_case1(infos, 0.0) + wx * ddx_f_plus_case1(infos, x_i) );
        infos->dydi_down[j_down*M + i] = infos->dxdi[i] * (1.0 + h_down) * ( w0 * ddx_f_minus_case1(infos, 0.0) + wx * ddx_f_minus_case1(infos, x_i) );

        infos->dydj_up[j_up*M + i]   = dhdj_up*( 0.5*H - ( w0 * f_plus_case1(infos, 0.0, infos->t) + wx * f_plus_case1(infos, x_i, infos->t) ) );
        infos->dydj_down[j_down*M + i] = dhdj_down*( 0.5*H + ( w0 * f_minus_case1(infos, 0.0, infos->t) + wx * f_minus_case1(infos, x_i, infos->t) ) );

        infos->D_up[j_up*M + i]   = infos->dxdi[i] * infos->dydj_up[j_up*M + i];
        infos->D_down[j_down*M + i] = infos->dxdi[i] * infos->dydj_down[j_down*M + i];

        infos->D_up_t_plus_1[j_up*M + i] = infos->dxdi[i]*dhdj_up*( 0.5*H - ( w0 * f_plus_case1(infos, 0.0, infos->t + infos->dt) + wx * f_plus_case1(infos, x_i, infos->t + infos->dt) ) );
        infos->D_down_t_plus_1[j_down*M + i] = infos->dxdi[i]*dhdj_down*( 0.5*H + ( w0 * f_minus_case1(infos, 0.0, infos->t + infos->dt) + wx * f_minus_case1(infos, x_i, infos->t + infos->dt) ) );
    }
}

void grow_grid(struct grid_and_data *infos){
    /*
        pre :
            infos : the structure we use
        post :
            grow the size of the grid, useful for supersonic cases
    */
    infos->L = 7.;
    infos->H = 6.*infos->c;
    //beta1 reste le même
    if (infos->resolution == 0.0){
        infos->beta_3 = 2.13188005; 
    }
    else if (infos->resolution == 0.5){
        infos->beta_3 = 2.11396953;
    }
    else if (infos->resolution == 1.0){
        infos->beta_3 = 2.10956048;
    }
    else{
       infos->beta_3 = 2.13188005;  
    }
}

struct grid_and_data init_grid_and_data(double M_inf, double grow_rate_M_inf, double t_end ,int number_of_wrotten_steps, double resolution, double amplitude_oscillation, double frequence_oscillation, double angle_of_attack,bool hdd_storage){
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
    infos->hdd_storage = hdd_storage;
    infos->resolution = resolution;


    if (resolution == 0.00){
        infos->M = 150; 
        infos->N = 30;  

        infos->beta_3 = 1.66109493; // déterminé va statement & python 
        infos->dt = 5e-6; 
    }
    else if (resolution == 0.5){
        infos->M = 450; 
        infos->N = 90;  

        infos->beta_3 = 1.64152192; // déterminé va statement & python  
        infos->dt = 1e-6;
    }
    else if (resolution == 1.){
        infos->M = 900; 
        infos->N = 360;  

        infos->beta_3 = 1.63669674; // déterminé va statement & python  
        infos->dt = 5e-7;
    }
    else{
        infos->M = 150; 
        infos->N = 30;  

        infos->beta_3 = 1.66109493; // déterminé va statement & python 
        infos->dt = 1e-5; 
    }
    infos->beta_1 = 1.08865949; // déterminé va statement & python, le même peu importe la résolution


    infos->M_inf_init = M_inf;
    infos->M_inf = M_inf;
    infos->grow_rate_M_inf = grow_rate_M_inf;

    // paramètres temporels
    infos->t = 0.;
    infos->t_end = t_end;
    infos->number_of_wrotten_steps = number_of_wrotten_steps;
    infos->write_interval_steps = t_end/(infos->dt*infos->number_of_wrotten_steps);

    if (M_inf <.6){
        infos->C = 0.;
        infos->eps_min = 1.;
    }
    else{
        infos->C = .1;
        infos->eps_min = 2.; 
    }


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
    infos->R = 3.*c;
    infos->H = 4.*c;
    infos->y_c = 0.0;
    infos->x_c = c/2.;

    infos->L = L;
    infos->L_u = 1.5;

    if (infos->M_inf > 1. || grow_rate_M_inf !=0.){
        grow_grid(infos);
    }

    infos->i_1 = (int)infos->M/6;
    infos->i_2 = (int)infos->M*5/6;
    infos->i_lead = 0;// déterminé  dans une autre fct
    infos->i_trail = 0; // déterminé  dans une autre fct

    infos->beta_plus = 1.; 
    infos->beta_minus = 1.;
    infos->beta = 1.;

    // paramètres de l'aile
    infos->AOA = angle_of_attack*M_PI/180.; // passe de deg à rad
    infos->AOA_current = 0.;
    infos->t_ramp = angle_of_attack*L/infos->c_sound_inf; //Proportionnel à l'ange d'attaque pour pas que ça aille trop vite et que ça pète
    
    if (angle_of_attack == 0.){
        infos->t_ramp = 67.67; // ça change rien peu importe t_ramp, AOA = 0 donc AOA_grow_rate = 0 donc on est bon simon, sinon ça crashe 
    }

    infos->AOA_grow_rate = 5.*infos->AOA/infos->t_end;

    infos->mu = amplitude_oscillation/2.;
    infos->omega = frequence_oscillation; 
    

    // grilles ,dérivées et caractéristiques
    infos->x = (double *)calloc(infos->M, sizeof(double));
    infos->y_up = (double *)calloc(infos->N*infos->M, sizeof(double));
    infos->y_down = (double *)calloc(infos->N*infos->M, sizeof(double));


    infos->dxdi = (double *)calloc(infos->M, sizeof(double));

    infos->dydi_up = (double *)calloc(infos->N*infos->M, sizeof(double));
    infos->dydj_up = (double *)calloc(infos->N*infos->M, sizeof(double));
    infos->dydt_up = (double *)calloc(infos->N*infos->M, sizeof(double));

    infos->dydi_down = (double *)calloc(infos->N*infos->M, sizeof(double));
    infos->dydj_down = (double *)calloc(infos->N*infos->M, sizeof(double));
    infos->dydt_down = (double *)calloc(infos->N*infos->M, sizeof(double));

    infos->D_up =(double *) calloc(infos->N*infos->M, sizeof(double));
    infos->D_down = (double *)calloc(infos->N*infos->M, sizeof(double));
    infos->D_up_t_plus_1 = (double *)calloc(infos->N*infos->M, sizeof(double));
    infos->D_down_t_plus_1 = (double *)calloc(infos->N*infos->M, sizeof(double));

    infos->s_up =  (double *)calloc(4*infos->N*infos->M, sizeof(double));
    infos->s_down =  (double *)calloc(4*infos->N*infos->M, sizeof(double));

    infos->s_star_up =  (double *)calloc(4*infos->N*infos->M, sizeof(double));
    infos->s_star_down =  (double *)calloc(4*infos->N*infos->M, sizeof(double));

    infos->F_i_up =  (double *)calloc(4*infos->N*infos->M, sizeof(double));
    infos->F_i_down =  (double *)calloc(4*infos->N*infos->M, sizeof(double));
    infos->F_i_star_up =  (double *)calloc(4*infos->N*infos->M, sizeof(double));
    infos->F_i_star_down =  (double *)calloc(4*infos->N*infos->M, sizeof(double));

    infos->F_j_up =  (double *)calloc(4*infos->N*infos->M, sizeof(double));
    infos->F_j_down =  (double *)calloc(4*infos->N*infos->M, sizeof(double));
    infos->F_j_star_up =  (double *)calloc(4*infos->N*infos->M, sizeof(double));
    infos->F_j_star_down =  (double *)calloc(4*infos->N*infos->M, sizeof(double));

    return *infos;
}

void free_data(struct grid_and_data *infos){
    /*
        pre :
            infos : the structure we use
        post :
           free all of the structure, we dont want it anymore ciao
    */
    free(infos->x);
    free(infos->dxdi);

    free(infos->y_up);
    free(infos->y_down);


    free(infos->D_up);
    free(infos->D_down);
    free(infos->D_up_t_plus_1);
    free(infos->D_down_t_plus_1);

    free(infos->dydi_up);
    free(infos->dydj_up);
    free(infos->dydt_up);

    free(infos->dydi_down);
    free(infos->dydj_down);
    free(infos->dydt_down);

    free(infos->s_up);
    free(infos->s_down);
    free(infos->s_star_up);
    free(infos->s_star_down);

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