import numpy as np
import matplotlib.pyplot as plt
import scipy as sc




#_____ betas finder

def beta_finder(beta_plus, L, L_u, c, M, i_1, i_2):
    y_c = 0
    def fct_beta_1(beta_1):
        return (L_u-c/2)*(1/i_1)*(2*beta_1/(np.sinh(2*beta_1))) - 2*c/(i_2-i_1)

    def fct_beta_3(beta_3):
        return (L-L_u-3*c/2)*(1/ (M-1 - i_2) )*(2*beta_3/(np.sinh(2*beta_3))) - 2*c/(i_2-i_1)

    def fct_beta_minus(beta_minus):
        test1 = (H/2 -y_c) *2*beta_minus/np.sinh(2*beta_minus) - (H/2 +y_c) *2*beta_plus/np.sinh(2*beta_plus)
        return beta_minus - beta_plus
    
    beta_1_guess = 0.5
    beta_1 = sc.optimize.fsolve(fct_beta_1, beta_1_guess)

    beta_3_guess = 0.5
    beta_3 = sc.optimize.fsolve(fct_beta_3, beta_3_guess)

    beta_minus = 1.
    beta_minus =  sc.optimize.fsolve(fct_beta_minus, beta_minus)

    

    return beta_1, beta_3, beta_minus

#____ grid 

def grid_x(beta_plus,L, L_u, c, M, i_1, i_2):
    
    beta_1, beta_3, _ = beta_finder(beta_plus,L, L_u, c, M, i_1, i_2)
    
    def g_upstream(i):
        g_i = (L_u - c/2)/L * np.tanh(beta_1*i/i_1)/np.tanh(beta_1)

        gprim_i_num = (L_u - c/2)*beta_1/i_1
        gprim_i_den = np.tanh(beta_1)*np.cosh(beta_1*i/i_1)*np.cosh(beta_1*i/i_1)

        return L*g_i-L_u, gprim_i_num/gprim_i_den
    
    def g_midstream(i):
        g_i = (L_u - c/2)/L + 2*c/L * (i-i_1)/(i_2 - i_1)
        gprim_i = 2*c/(i_2-i_1)
        return L*g_i -L_u, gprim_i
    
    def g_downstream(i):
        term1 = (L_u + 3*c/2)/L
        
        term2 = 1 - term1
        arg = beta_3* (1- (i - i_2)/(M-1-i_2))
        term3 = 1 - ( np.tanh ( arg)/np.tanh(beta_3))

        g_i = term1 + term2*term3
    
        f_i = 1-(i-i_2)/(M-1-i_2)
        fprim_i = -1/(M-1-i_2)
        gprim_i_num = -L*(term2)*beta_3*fprim_i


        gprim_i_den = np.tanh(beta_3)*np.cosh(beta_3*f_i)*np.cosh(beta_3*f_i)
        return L*g_i-L_u, gprim_i_num/gprim_i_den
    
    x = np.zeros(M)
    dxdi = np.zeros(M)


    for i in range(0, i_1):
        x[i], dxdi[i] = g_upstream(i)

    for i in range(i_1, i_2):
        x[i], dxdi[i] = g_midstream(i)

    for i in range(i_2, M):
        x[i], dxdi[i] = g_downstream(i)

    

    return x, dxdi

def grid_y(beta_plus,beta_minus,x ,M, N, c, R,H,  mu, omega, t):
    AoA = 5*np.pi/180
    def f_plus(x_i):
        y = np.sqrt(-(x_i-x_c)**2+R**2)-np.sqrt(-x_c**2+R**2) + (1-2*x_i/c)*  mu*c*np.sin(omega*t) +  (1-2*x_i/c)*c*np.tan(AoA)
        return y
    
    def f_minus(x_i):
        y = -np.sqrt(-(x_i-x_c)**2+R**2)+np.sqrt(-x_c**2+R**2) + (1-2*x_i/c)*  mu*c*np.sin(omega*t) +  (1-2*x_i/c)*c*np.tan(AoA)
        return y
    
    def circle_bounds_finder(x):
        index_0 = 0
        index_1 = 0
        while x[0][index_0] <= 1e-6:
            index_0 +=1

        index_1 = index_0

        while (x[0][index_1]) <=1-1e-6:
            index_1 +=1

        return index_0, index_1

    beta = 1
    y_c = 0
    x_c = c/2

    y_up = np.zeros((N,M))
    y_down = np.zeros((N,M))
    index_x_0, index_x_1 = circle_bounds_finder(x)


    # grid iniriale
    for j in range(N):
        arg_up =  1 - j/(N-1)
        arg_down =   j/(N-1)

        hj_up = 1 - np.tanh(beta_plus*arg_up)/np.tanh(beta)
        hj_down = np.tanh(beta_minus*arg_down)/np.tanh(beta)-1

        for i in range(0, index_x_0):
            y_up[j][i]=     H/2*hj_up + arg_up*f_plus(0)
            y_down[j][i] =   H/2*hj_down + arg_down*f_minus(0)

        for i in range(index_x_0, index_x_1):
            y_up[j][i] =    H/2*hj_up + arg_up*f_plus(x[j][i])
            y_down[j][i] =  H/2*hj_down + arg_down*f_minus(x[j][i])

        for i in range(index_x_1, M):
            y_up[j][i]=     H/2*hj_up + arg_up*f_plus(1)
            y_down[j][i] =   H/2*hj_down + arg_down*f_minus(1)

    return y_up,y_down

def grid_plotter( beta_plus,L, L_u, c,R, H,  M, N,i_1, i_2,  t, mu, omega):

    _,_, beta_minus= beta_finder(beta_plus, L, L_u, c, M, i_1, i_2)
    x_half, _ = grid_x(beta_plus,L, L_u, c, M, i_1, i_2)
    x_half = np.repeat(x_half, N).reshape(M,N).T
    y_up, y_down = grid_y(beta_plus,beta_minus,x_half, M, N, c, R, H, mu, omega, t)

    plt.pcolormesh(x_half, y_up, np.zeros((N-1, M-1)), facecolor='none', edgecolors='red', linewidth=0.5)
    plt.pcolormesh(x_half, y_down,np.zeros((N-1, M-1)), facecolor='none', edgecolors='black', linewidth=0.5)

    plt.show()

#_____ CFL

def get_max_dt_using_cfl_condition_in_2D(beta_plus, L, L_u, c,R, H,  M, N,i_1, i_2,  t, mu, omega, max_speed):

    _,_, beta_minus= beta_finder(beta_plus, L, L_u, c, M, i_1, i_2)
     
    x_half, _ = grid_x(beta_plus,L, L_u, c, M, i_1, i_2)
    x_half = np.repeat(x_half, N).reshape(M,N).T
    y_up, _ = grid_y(beta_plus,beta_minus,x_half, M, N, c, R, H, mu, omega, t)

    dx = np.min(np.diff(x_half, axis=1))
    dy = np.min(np.diff(y_up, axis=0))

    dt_cfl_x = dx / max_speed
    dt_cfl_y = dy / max_speed

    dt_cfl = min(dt_cfl_x, dt_cfl_y)

    return dt_cfl

#_____ fixed parameters

c = 1
R = 3*c
H = 6*c

L = 7
L_u = 1.5

# low resolution 

M = 150
N = 30


# mid resol

M = 450
N = 90

# high parce que pq pas


M = 900
N = 360

# ouais
beta_plus = 1
i_1 = M//6
t = 1.0
i_2 = 5*M//6
mu = 0.1


omega = np.pi*0.5/t
omega = 0

mach = 1.6
max_speed = mach*340 #m/s vaut mach 5

beta_1, beta_3, beta_minus = beta_finder(beta_plus,L, L_u, c, M, i_1, i_2)

print("beta_1: ", beta_1)
print("beta_3: ", beta_3)


grid_plotter(beta_plus,L, L_u, c, R, H, M, N, i_1, i_2,t, mu, omega)


#cfl_dt = get_max_dt_using_cfl_condition_in_2D(beta_plus,L, L_u, c,R, H,  M, N,i_1, i_2,  t, mu, omega, max_speed)
#print(cfl_dt)