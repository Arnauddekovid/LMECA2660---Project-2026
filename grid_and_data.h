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
    double M_inf;
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
    double* angle;
    double y_c;
    double x_c;

    double airfoil; // True si on teste avec l'aile, false sinon
    bool print_evolution_of_mass;
    bool debug_solver; // enable detailed debug prints (useful for failures near airfoil)

    double L;
    double L_u;
    //double L_d;

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

    // grids
    double* x;  
    double**y_up;
    double**y_up_t_plus_1; // the grid at the t + dt in order to get the tempo dervative
    double**y_down;
    double**y_down_t_plus_1; 

    double ** mac_cormack_f_up_checker;
    double ** mac_cormack_f_down_checker;

    double ** mac_cormack_b_up_checker;
    double ** mac_cormack_b_down_checker;


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

void leading_and_trailing_edge_finder(struct grid_and_data *infos){
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

double f_plus(struct grid_and_data*infos, int index){   
    double x_i = infos->x[index];
    double value = sqrt(-(x_i-infos->x_c)*(x_i-infos->x_c)+infos->R*infos->R)- sqrt(-infos->x_c*infos->x_c+infos->R*infos->R);

    return value /2.;
}

double f_minus(struct grid_and_data*infos, int index){
    double x_i = infos->x[index];
    double value = sqrt(-(x_i-infos->x_c)*(x_i-infos->x_c)+infos->R*infos->R)- sqrt(-infos->x_c*infos->x_c+infos->R*infos->R);
    return - value /2.;
}

double f_plus_fs(struct grid_and_data*infos, int index){
    double x_i = infos->x[index];
    return (1. - cos(2*M_PI*x_i))/100.;
}

double f_minus_fs(struct grid_and_data*infos, int index){
    double x_i = infos->x[index];
    return -(1. - cos(2*M_PI*x_i))/100.;
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
        infos->x[i] = ( (infos->L_u - infos->c/2) + 2*infos->c * (i-infos->i_1)/(infos->i_2 - infos->i_1))- infos->L_u;
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

    leading_and_trailing_edge_finder(infos);   
    
    for (int j =0; j <infos->N;j++){
        for (int i = 0; i < infos->M; i++){
            infos->y_up_t_plus_1[j][i]   = (infos->H/2)*(  1 - tanh(infos->beta_plus * ( 1. - (double)j/(infos->N-1)))/tanh(infos->beta));
            infos->y_down_t_plus_1[j][i] = (infos->H/2)*(-(1 - tanh(infos->beta_minus* ((double)j/(infos->N-1)))/tanh(infos->beta)));
        }
        for (int i = infos->i_lead; i < infos->i_trail+1; i++){
            infos->y_up_t_plus_1[j][i] += infos->airfoil*(infos->H/2)*(1.-(double) j/(infos->N-1))*f_plus(infos, i);
            infos->y_down_t_plus_1[j][i] += infos->airfoil*(infos->H/2)*((double)j/(infos->N-1))*f_minus(infos, i);



        }
    }
}

void update_y_grid(struct grid_and_data *infos){

    // nous fait la grid selon y et dy/dmachin

    int i_lead = infos->i_lead;
    int i_trail = infos->i_trail;

    int di = 2;

    // on fait un smoothing entre i_lead -di jusqu'a i_lead + di soit 2*di+1 points ouy on smooth pour pas que les dérivées partent en steak

    for (int j =0; j <infos->N;j++){

        double h_j_up   = 1. - tanh(infos->beta_plus*(1.-(double)j/ ((double)(infos->N-1))))/tanh(infos->beta);
        double h_j_down = tanh(infos->beta_minus*((double)j/ ((double)(infos->N-1))))/tanh(infos->beta)-1.;
        
        double numerator_up = infos->beta_plus/((double)(infos->N-1));
        double denominator_up = tanh(infos->beta)*cosh(infos->beta_plus * ( 1. - (double)j/(infos->N-1)))*cosh(infos->beta_plus * ( 1. - (double)j/(infos->N-1)));

        double dh_updj_coef = (infos->H/2)*numerator_up/denominator_up;

        double numerator_down = infos->beta_minus/((double)(infos->N-1));
        double denominator_down = tanh(infos->beta)*cosh(infos->beta_minus * ( (double)j/(infos->N-1)))*cosh(infos->beta_minus * ( (double)j/(infos->N-1)));

        double dh_downdj_coef = (infos->H/2)*numerator_down/denominator_down;

        for (int i = 0; i < infos->M; i++){
            infos->y_up[j][i]   = infos->y_up_t_plus_1[j][i];
            infos->y_down[j][i] = infos->y_down_t_plus_1[j][i];

            infos->y_up_t_plus_1[j][i]   = (infos->H/2)*(  1. - tanh(infos->beta_plus * ( 1. - (double)j/(infos->N-1.)))/tanh(infos->beta));
            infos->y_down_t_plus_1[j][i] = (infos->H/2)*(-(1. - tanh(infos->beta_minus* ( (double)j/(infos->N-1.)))/tanh(infos->beta)));
            
            infos->dy_updj[j][i] = dh_updj_coef;
            infos->dy_downdj[j][i] = dh_downdj_coef;

            infos->dy_updi[j][i] = 0.;
            infos->dy_downdi[j][i] = 0.;
        }

        //leading edge smoothing

        for (int i = infos->i_lead;i< infos->i_lead+di;i++){
            double step = (double)(i - infos->i_lead) / (double)di;
            double scale = step*step*(3.-2.*step);


            infos->y_up_t_plus_1[j][i] += scale*infos->airfoil*(infos->H/2)*(1.-(double) j/(infos->N-1))*f_plus(infos, i);
            infos->y_down_t_plus_1[j][i] += scale*infos->airfoil*(infos->H/2)*(1. - (double) (infos->N-1-j)/(infos->N-1))*f_minus(infos, i);

            infos->dy_updj[j][i]   += scale*infos->airfoil*dh_updj_coef*f_plus(infos, i);
            infos->dy_downdj[j][i] -= scale*infos->airfoil*dh_downdj_coef*f_minus(infos, i);

            double df_plusdi = -0.5*((infos->x[i]-infos->x_c)/sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c)));
            double df_minusdi = -df_plusdi;
        
            infos->dy_updi[j][i] += scale*infos->airfoil*(infos->H/2)*infos->dxdi[i]*(1.-h_j_up)*df_plusdi;
            infos->dy_downdi[j][i] += scale*infos->airfoil*(infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }


        for (int i = infos->i_lead+di; i <= infos->i_trail-di; i++){
        
            infos->y_up_t_plus_1[j][i] += infos->airfoil*(infos->H/2)*(1.-(double) j/(infos->N-1))*f_plus(infos, i);
            infos->y_down_t_plus_1[j][i] += infos->airfoil*(infos->H/2)*(1. - (double) (infos->N-1-j)/(infos->N-1))*f_minus(infos, i);

            infos->dy_updj[j][i]   += infos->airfoil*dh_updj_coef*f_plus(infos, i);
            infos->dy_downdj[j][i] -= infos->airfoil*dh_downdj_coef*f_minus(infos, i);

            double df_plusdi = -0.5*((infos->x[i]-infos->x_c)/sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c)));
            double df_minusdi = -df_plusdi;
        
            infos->dy_updi[j][i] += infos->airfoil*(infos->H/2)*infos->dxdi[i]*(1.-h_j_up)*df_plusdi;
            infos->dy_downdi[j][i] += infos->airfoil*(infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;

        }
        //in the beginning of the airfoil we add


        for (int i = infos->i_trail-di;i< infos->i_trail;i++){
            double step = (double)(infos->i_trail - i) / (double)di; 
            double scale = step*step*(3.-2.*step);
            infos->y_up_t_plus_1[j][i] += scale*infos->airfoil*(infos->H/2)*(1.-(double) j/(infos->N-1))*f_plus(infos, i);
            infos->y_down_t_plus_1[j][i] += scale*infos->airfoil*(infos->H/2)*(1. - (double) (infos->N-1-j)/(infos->N-1))*f_minus(infos, i);

            infos->dy_updj[j][i]   += scale*infos->airfoil*dh_updj_coef*f_plus(infos, i);
            infos->dy_downdj[j][i] -= scale*infos->airfoil*dh_downdj_coef*f_minus(infos, i);

            double df_plusdi = -0.5*((infos->x[i]-infos->x_c)/sqrt(infos->R*infos->R - (infos->x[i]-infos->x_c)*(infos->x[i]-infos->x_c)));
            double df_minusdi = -df_plusdi;
        
            infos->dy_updi[j][i] += scale*infos->airfoil*(infos->H/2)*infos->dxdi[i]*(1.-h_j_up)*df_plusdi;
            infos->dy_downdi[j][i] += scale*infos->airfoil*(infos->H/2)*infos->dxdi[i]*(1.+h_j_down)*df_minusdi;
        }
        


        for (int i = 0; i < infos->M; i++){
            infos->dy_updt[j][i] = (infos->y_up_t_plus_1[j][i] - infos->y_up[j][i] )/infos->dt;
            infos->dy_downdt[j][i] = (infos->y_down_t_plus_1[j][i] - infos->y_down[j][i])/infos->dt;

            infos->D_up[j][i] = infos->dxdi[i]*infos->dy_updj[j][i];
            infos->D_down[j][i] = infos->dxdi[i]*infos->dy_downdj[j][i];
        }
    }
}

struct grid_and_data init_grid_and_data(double M_inf, double t_end, double dt,int step_max, double C, double eps_min){

    struct grid_and_data infos_local;

    struct grid_and_data *infos = &infos_local;

    infos->airfoil = true; // true si on prend en compte l'airfoil pour les calculs

    infos->M_inf = M_inf;
    infos->C = C;
    infos->eps_min = eps_min;

    infos->gamma = 1.4;
    infos->T_inf = 300.0;
    infos->p_inf = 1e5;
    infos->R_fluid = 287.1;

    infos->rho_inf =infos->p_inf/ (infos->R_fluid *infos->T_inf);
    infos->c_sound_inf = sqrt(infos->gamma*infos->R_fluid*infos->T_inf);

    infos->T_0_inf = infos->T_inf*(1.+0.5*(infos->gamma-1.)*infos->M_inf*infos->M_inf);
    infos->p_0_inf = infos->p_inf*pow(infos->T_0_inf / infos->T_inf, infos->gamma/(infos->gamma-1) );

    double c = 1.;
    double L = 5.;

    infos->c = c;
    infos->R = 3*c;
    infos->H = 4*c;
    infos->y_c = 0.0;
    infos->x_c = c/2;

    infos->L = L;
    infos->L_u = 1.5;
    //infos->L_d = L - (infos->L_u-c);

    infos->M = 150;// 150 au début
    infos->i_1 = (int)infos->M/5;
    infos->i_2 = (int)infos->M*4/5;
    infos->i_lead = 0;
    infos->i_trail = 0;
    infos->N = 30; //30 au débur

    infos->t = 0.;
    infos->dt = dt; 
    infos->t_end = t_end;
    infos->step_max = step_max;

    infos->beta_1 = 0.81106561;
    infos->beta_3 = 1.44500835;

    infos->beta_plus = 1.;
    infos->beta_minus = 1.;
    infos->beta = 1.;
    
    infos->angle = (double *)calloc( (int)(t_end/infos->dt), sizeof(double));
    infos->x = (double *)calloc(infos->M, sizeof(double));

    infos->y_up = (double **)calloc(infos->N, sizeof(double*));
    infos->y_down = (double **)calloc(infos->N, sizeof(double*));

    infos->y_up_t_plus_1 = (double **)calloc(infos->N, sizeof(double*));
    infos->y_down_t_plus_1 = (double **)calloc(infos->N, sizeof(double*));

    infos->mac_cormack_f_up_checker = (double **)calloc(infos->N, sizeof(double*));
    infos->mac_cormack_f_down_checker = (double **)calloc(infos->N, sizeof(double*));

    infos->mac_cormack_b_up_checker = (double **)calloc(infos->N, sizeof(double*));
    infos->mac_cormack_b_down_checker = (double **)calloc(infos->N, sizeof(double*));


    infos->dxdi = (double *)calloc(infos->M, sizeof(double));

    infos->dy_updi = (double **)calloc(infos->N, sizeof(double*));
    infos->dy_updj = (double **)calloc(infos->N, sizeof(double*));
    infos->dy_updt = (double **)calloc(infos->N, sizeof(double*));

    infos->dy_downdi = (double **)calloc(infos->N, sizeof(double*));
    infos->dy_downdj = (double **)calloc(infos->N, sizeof(double*));
    infos->dy_downdt = (double **)calloc(infos->N, sizeof(double*));

    infos->D_up = (double **)calloc(infos->N, sizeof(double*));
    infos->D_down = (double **)calloc(infos->N, sizeof(double*));


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

        infos->mac_cormack_f_up_checker[j] = (double*)calloc(infos->M, sizeof(double));
        infos->mac_cormack_f_down_checker[j] = (double*)calloc(infos->M, sizeof(double));

        infos->mac_cormack_b_up_checker[j] = (double*)calloc(infos->M, sizeof(double));
        infos->mac_cormack_b_down_checker[j] = (double*)calloc(infos->M, sizeof(double));

        infos->D_up[j] = (double*)calloc(infos->M, sizeof(double));
        infos->D_down[j] = (double*)calloc(infos->M, sizeof(double));
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

        free(infos->mac_cormack_f_up_checker[j]);
        free(infos->mac_cormack_f_down_checker[j]);

        free(infos->mac_cormack_b_up_checker[j]);
        free(infos->mac_cormack_b_down_checker[j]);

        free(infos->D_up[j]);
        free(infos->D_down[j]);

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

    free(infos->mac_cormack_f_up_checker);
    free(infos->mac_cormack_f_down_checker);

    free(infos->mac_cormack_b_up_checker);
    free(infos->mac_cormack_b_down_checker);

    free(infos->D_up);
    free(infos->D_down);

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