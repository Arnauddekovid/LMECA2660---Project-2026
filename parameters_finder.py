import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

#_____ fixed parameters

c = 1
R = 3*c
H = 4*c

#_____ our parameters 
L = 5
L_u = 1.5
L_d = L - (L_u-c)

M = 150
i_1 = M//5
i_2 = 4*M//5

N = 30


# ,new parametres

M = 450
i_1 = M//5
i_2 = 4*M//5

N = 90



mu = 0.5
omega = np.pi/2*1/0.1

airfoil = True

#_____ betas finder

def fct_beta_1(beta_1):
    return (L_u-c/2)*(1/i_1)*(2*beta_1/(np.sinh(2*beta_1))) - 2*c/(i_2-i_1)

def fct_beta_3(beta_3):
    return (L-L_u-3*c/2)*(1/ (M-1 - i_2) )*(2*beta_3/(np.sinh(2*beta_3))) - 2*c/(i_2-i_1)


beta_1_guess = 0.5
beta_1 = sc.optimize.fsolve(fct_beta_1, beta_1_guess)

beta_3_guess = 0.5
beta_3 = sc.optimize.fsolve(fct_beta_3, beta_3_guess)

#____ functions 

def grid_x():
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

def grid_y_before(x,t):

    def f_plus(x_i,t):
        y = np.sqrt(-(x_i-x_c)**2+R**2)-np.sqrt(-x_c**2+R**2) + (1-2*x_i/c)*mu*c*np.sin(omega*t)
        return y*0.5
    
    def f_minus(x_i,t):
        y = -np.sqrt(-(x_i-x_c)**2+R**2)+np.sqrt(-x_c**2+R**2) + (1-2*x_i/c)*mu*c*np.sin(omega*t)
        return y*0.5
    
    def circle_bounds_finder(x):
        index_0 = 0
        index_1 = 0
        while x[0][index_0] <=0:
            index_0 +=1

        index_1 = index_0

        while x[0][index_1] <=1:
            index_1 +=1

        return index_0, index_1


    beta_minus= 1
    beta_plus = 1
    beta = 1
    y_c = 0
    x_c = c/2

    y_up = np.zeros((N,M))
    y_down = np.zeros((N,M))
    index_x_0, index_x_1 = circle_bounds_finder(x)

    time_coef = mu*c*np.sin(omega*t)

    # grid iniriale
    for j in range(N):
        arg_up =  1 - j/(N-1)
        arg_down =   j/(N-1)

        coef_up = 1 - np.tanh(beta_plus*arg_up)/np.tanh(beta)
        coef_down = 1 - np.tanh(beta_minus*arg_down)/np.tanh(beta)

        for i in range(0, index_x_0):
            y_up[j][i]=     (H/2- time_coef)*coef_up + time_coef
            y_down[j][i] = -(H/2 +time_coef)*coef_down +time_coef

        for i in range(index_x_0, index_x_1+1):
            y_up[j][i] = H/2*(coef_up + arg_up*f_plus(x[j][i], t))
            y_down[j][i] =  -H/2*(coef_down - arg_down*f_minus(x[j][i], t))

        for i in range(index_x_1+1, M):
            y_up[j][i]= (H/2 +time_coef)*coef_up  - time_coef
            y_down[j][i] = -(H/2-time_coef)*coef_down -time_coef

    
    
    return y_up,y_down

def get_max_dt_using_cfl_condition_in_2D(x, y, max_speed):
    dx = np.min(np.diff(x, axis=1))
    dy = np.min(np.diff(y, axis=0))

    dt_cfl_x = dx / max_speed
    dt_cfl_y = dy / max_speed

    dt_cfl = min(dt_cfl_x, dt_cfl_y)

    return dt_cfl


decoupled = False

x_half, dxdi_half = grid_x()

print("beta_1:", beta_1)
print("beta_3:", beta_3)


"""
for i in range(len(dxdi_half)-1):
    print(i, float(dxdi_half[i]), float(x_half[i+1]-x_half[i]))
"""
t = 0.00 #s
x_half = np.repeat(x_half, N).reshape(M,N).T
y_up, y_down = grid_y_before(x_half, t)

#plt.scatter(x, y, s=1)
#plt.pcolormesh(x_half, y_up, np.zeros((N-1, M-1)), facecolor='none', edgecolors='red', linewidth=0.5)
plt.pcolormesh(x_half, y_down,np.zeros((N-1, M-1)), facecolor='none', edgecolors='black', linewidth=0.5)

plt.show()

cfl_dt = get_max_dt_using_cfl_condition_in_2D(x_half, y_up, 1000)

print(cfl_dt)