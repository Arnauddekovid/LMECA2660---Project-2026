import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os 
from scipy.ndimage import gaussian_filter  

def convert_time_to_readable_format():
    """
    pre :
        none
    post:
        returns the current time in a human readable format (day month year hour:minute:second)
    """

    time_struct = time.localtime()
    readable_time = time.strftime("%d %B %Y %H:%M:%S", time_struct)
    return readable_time

def extract_data(filename):
    dico_data ={}

    x = []
    t_liste = []
    M_inf_liste = []

    y_all = []
    s_0_all = []
    s_1_all = []
    s_2_all = []
    s_3_all = []

    s_all = []


    with open(filename, 'r') as f:

        file = f.readlines()
        infos_str = file[0].split(" ")[:-1]

        number_of_data = len(infos_str)//3

        for i in range(number_of_data):
            var = infos_str[3*i]
            value = infos_str[3*i+2]
            dico_data[var] = float(value)

        x_str = file[1].split(" ")[:-1]

        t = float(file[2].split(" ")[4]),
        M_inf = float(file[2].split(" ")[7])

        for x_j in x_str:
            x.append(float(x_j))
        i = 3
        np.array(x)
        M = 0
        N = 0

        x_reshape_already_done = False


        while file[i-1].split(" ")[0] != "end" :


            y = []
            s_0 = []
            s_1 = []
            s_2 = []
            s_3 = []
            
            ligne_i =file[i-1].split(" ")
            t_i = ligne_i[4]
            M_inf_i =ligne_i[7]

            t_liste.append(float(t_i))
            M_inf_liste.append(float(M_inf_i))

            # y d'abord

            while (len(file[i].split(" "))> 10):
                y_i = file[i].split(" ")
                y_liste = []
                for y_ij in y_i[:-1]:
                    y_liste.append(float(y_ij))
                y.append(np.array(y_liste))
                i+=1
            i+=1

            # s_0

            while (len(file[i].split(" "))> 10):
                u_i = file[i].split(" ")
                u_liste = []
                for u_ij in u_i[:-1]:

                    u_liste.append(float(u_ij))
                s_0.append(np.array(u_liste))
                i+=1 
            i+=1

            #s_1
            while (len(file[i].split(" "))>10):
                u_i = file[i].split(" ")
                u_liste = []
                for u_ij in u_i[:-1]:

                    u_liste.append(float(u_ij))
                s_1.append(np.array(u_liste))
                i+=1 

            i+=1
            #s_2

            while (len(file[i].split(" "))>10):
                u_i = file[i].split(" ")
                u_liste = []
                for u_ij in u_i[:-1]:

                    u_liste.append(float(u_ij))
                s_2.append(np.array(u_liste))
                i+=1 
            i+=1
            #s_3

            while (len(file[i].split(" "))>10):
                u_i = file[i].split(" ")
                u_liste = []
                for u_ij in u_i[:-1]:

                    u_liste.append(float(u_ij))
                s_3.append(np.array(u_liste))
                i+=1 
            i+=1


            np.array(y)
            np.array(s_0)
            np.array(s_1)
            np.array(s_2)
            np.array(s_3)
 

            if not x_reshape_already_done:
                M = len(x)
                N = len(y)
                x = np.repeat(x, len(y)).reshape(M,N).T
                x_reshape_already_done = True
            
            y = np.reshape(y, (N, M))
            s_0 = np.reshape(s_0, (N, M))   
            s_1 = np.reshape(s_1, (N, M))   
            s_2 = np.reshape(s_2, (N, M))   
            s_3 = np.reshape(s_3, (N, M))   

            y_all.append(y)
            s_0_all.append(s_0)
            s_1_all.append(s_1)
            s_2_all.append(s_2)
            s_3_all.append(s_3)


    s_all = np.array([np.array(s_0_all), np.array(s_1_all),np.array(s_2_all), np.array(s_3_all)])

    
    return dico_data,  np.array(t_liste),M_inf_liste, x,y_all, s_all

def evolution_animation_grid(filename):

    _, t_liste,_, x, y_all, _= extract_data(filename)
    N = len(x)//2
    M = len(x[0])

    fig, ax = plt.subplots()

    mesh_up = ax.pcolormesh(
        x[N:],
        y_all[0][N:],
       np.zeros((N-1,M-1)),facecolor='none', edgecolors='red', linewidth=0.2)

    cbar = plt.colorbar(mesh_up, ax=ax)
    title = ax.set_title(f"t = {t_liste[0]:.3f}")

    def update(frame):
        ax.clear()
        mesh_up = ax.pcolormesh(x[N:],y_all[frame][N:],np.zeros((N-1,M-1)),facecolor='none', edgecolors='red', linewidth=0.2)
        mesh_down = ax.pcolormesh(x[N:],y_all[frame][:N],np.zeros((N-1,M-1)),facecolor='none', edgecolors='black', linewidth=0.2)
        ax.set_title(filename[:11]+f" t = {t_liste[frame]:.6f}")
        ax.axis('equal')

        return mesh_up,mesh_down

    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(t_liste),
        interval=100,
        blit=False
    )

    
    time = convert_time_to_readable_format() 

    string = "Animations/grid "+filename[:-4] + " at " + time + ".mp4"

    anim.save(string, writer='ffmpeg', fps=30)  

#_________

def data_over_time_animation(datas, data_name, plot_title,  x, y_all, M_inf_liste, t_liste):
    vmin = min(np.min(data) for data in datas)
    vmax = max(np.max(data) for data in datas)
    N = len(x)//2
    fig, ax = plt.subplots()

    mesh_up = ax.pcolormesh(
        x[N:],
        y_all[0][N:],
        datas[0][N:],
        vmin=vmin,
        vmax=vmax,            
        shading = "gouraud"
    )

    mesh_down = ax.pcolormesh(
        x[N:],
        y_all[0][:N],
        datas[0][:N],
        vmin=vmin,
        vmax=vmax,
        shading = "gouraud"

    )

    cbar = plt.colorbar(mesh_up, ax=ax)
    title = ax.set_title(f"{plot_title} with M_inf = {M_inf_liste[0]:.3f} t = {t_liste[0]:.3f}")
    ax.axis('equal')

    # --- Update function ---
    def update(frame):
        ax.clear()

        mesh_up = ax.pcolormesh(
            x[N:],
            y_all[frame][N:],
            datas[frame][N:],
            vmin=vmin,
            vmax=vmax,
            shading = "gouraud"
            
        )

        mesh_down = ax.pcolormesh(
            x[N:],
            y_all[frame][:N],
            datas[frame][:N],
            vmin=vmin,
            vmax=vmax,
            shading = "gouraud"

        )

        ax.set_ylim(-2, 2)

        # Puis impose le ratio égal sans toucher aux limites
        ax.set_aspect('equal', adjustable='box')

        ax.set_title(f"{plot_title} with M_inf = {M_inf_liste[0]:.3f} t = {t_liste[0]:.3f}")

        return mesh_up,mesh_down

    # --- Animation ---
    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(t_liste),
        interval=100,
        blit=False
    )

    # --- Save ---
    
    time = convert_time_to_readable_format() 

    #string = "Animations/wc over U_inf with "+filename[:-4].replace("_", " ") + " at " + time + ".mp4"
    string = "Animations/" + data_name+  " at " + time + ".mp4"

    anim.save(string, writer='ffmpeg', fps=30)  

def s_over_time_animation(filename):

    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename)

    rho_all = s_all[0]
    u_all = s_all[1]/rho_all
    v_all = s_all[2]/rho_all
    U0_all = s_all[3]/rho_all

    trucs = [rho_all, u_all, v_all, U0_all]

    for i in range(4):

        s_i = trucs[i]
        title_plot = "s" + str(i)
        data_name = title_plot + " "+  filename[:-4].replace("_", " ")
        data_over_time_animation(s_i, data_name, title_plot, x, y_all, M_inf_liste, t_liste)

def s0_over_time_animation(filename):
    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename)
    rho_all = s_all[0]
    title_plot = "rho_all" 
    data_name = title_plot + " "+  filename[:-4].replace("_", " ")
    data_over_time_animation(rho_all, data_name, title_plot, x, y_all, M_inf_liste, t_liste)
    
def Mach_number_over_time_animation(filename):

    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename)

    gamma = dico_data["gamma"]
    rho_all = s_all[0]
    u_all = s_all[1]/rho_all
    v_all = s_all[2]/rho_all
    U0_all = s_all[3]/rho_all

    p_all_over_rho_all = (gamma-1)*(U0_all -0.5*(u_all*u_all + v_all*v_all))

    Mach_all = np.sqrt( (u_all*u_all + v_all*v_all)/(gamma*p_all_over_rho_all))


    #-----
    title_plot = "Mach number"
    data_name = title_plot + " "+  filename[:-4].replace("_", " ")
    data_over_time_animation(Mach_all, data_name, title_plot, x, y_all, M_inf_liste, t_liste)

def Cp_over_time_animation(filename):
    # num : p_over_p_inf -1
    # den : gamma Mach_inf**2/2
    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename)

    gamma = dico_data["gamma"]
    M_inf = dico_data["M_inf"]
    p_0_inf = dico_data["p_0_inf"]
    T_0_inf = dico_data["T_0_inf"]

    rho_all = s_all[0]
    u_all = s_all[1]/rho_all
    v_all = s_all[2]/rho_all
    U0_all = s_all[3]/rho_all

    T_inf =  T_0_inf/(1+(gamma-1)*M_inf*M_inf/2)
    p_inf = p_0_inf*(T_inf/T_0_inf)**(gamma/(gamma-1))


    p_all = (gamma-1)*rho_all*(U0_all -0.5*(u_all*u_all + v_all*v_all))
    
    num = p_all/p_inf-1
    den = gamma*M_inf*M_inf/2
    Cp_all = num/den

    #-----
    title_plot = "Cp"
    data_name = "Cp "+ filename[:-4].replace("_", " ")
    data_over_time_animation(Cp_all, data_name, title_plot, x, y_all, M_inf_liste, t_liste)

def normalized_total_enthalpy_over_time_animation(filename):
    # H_0/H_0_inf
    # H_0 = U_0 +p/rho = p/[(gamma-1)rho] +0.5(u**2*v**2) +p/rho
    #                  = gamma/(gamma-1)*p/rho  +0.5(u**2*v**2)
    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename)

    gamma = dico_data["gamma"]
    M_inf = dico_data["M_inf"]
    R_fluid = dico_data["R_fluid"]
    p_0_inf = dico_data["p_0_inf"]
    T_0_inf = dico_data["T_0_inf"]

    T_inf =  T_0_inf/(1+(gamma-1)*M_inf*M_inf/2)
    p_inf = p_0_inf*(T_inf/T_0_inf)**(gamma/(gamma-1))
    rho_inf = p_inf/(R_fluid*T_inf)
    u_inf = M_inf*np.sqrt(gamma*R_fluid*T_inf)

    H_0_inf = gamma/(gamma-1) * p_inf/rho_inf+0.5*u_inf*u_inf

    rho_all = s_all[0]
    u_all = s_all[1]/rho_all
    v_all = s_all[2]/rho_all
    U0_all = s_all[3]/rho_all
    p_all_over_rho_all = (gamma-1)*(U0_all - 0.5*(u_all*u_all + v_all*v_all))

    H_0_all =  gamma/(gamma-1)*p_all_over_rho_all+ 0.5*(u_all*u_all + v_all*v_all)

    H_0_over_H_0_inf = H_0_all/H_0_inf

    #-----
    title_plot = "H_0/H_0_inf"
    data_name = "Total enthalpy "+filename[:-4].replace("_", " ")
    data_over_time_animation(H_0_over_H_0_inf, data_name, title_plot, x, y_all, M_inf_liste, t_liste)

def normalized_total_pressure_over_time_animation(filename):
    # p_0/p_0_inf
    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename)

    gamma = dico_data["gamma"]
    M_inf = dico_data["M_inf"]
    R_fluid = dico_data["R_fluid"]
    p_0_inf = dico_data["p_0_inf"]
    T_0_inf = dico_data["T_0_inf"]

    T_inf =  T_0_inf/(1+(gamma-1)*M_inf*M_inf/2)
    p_inf = p_0_inf*(T_inf/T_0_inf)**(gamma/(gamma-1))
    rho_inf = p_inf/(R_fluid*T_inf)
    u_inf = M_inf*np.sqrt(gamma*R_fluid*T_inf)

    rho_all = s_all[0]
    u_all = s_all[1]/rho_all
    v_all = s_all[2]/rho_all
    U0_all = s_all[3]/rho_all
    p_all = (gamma-1)*rho_all*(U0_all - 0.5*(u_all*u_all + v_all*v_all))
    Mach_all = np.sqrt( (u_all*u_all + v_all*v_all)/(gamma*p_all/rho_all))

    p_0_all = p_all*(1+(gamma-1)*Mach_all*Mach_all/2)**(gamma/(gamma-1))

    p_O_over_p_0_inf = p_0_all/p_0_inf

    title_plot = "p_0/p_0_inf"
    data_name =  "Total pressure " +  filename[:-4].replace("_", " ")
    data_over_time_animation(p_O_over_p_0_inf, data_name, title_plot, x, y_all, M_inf_liste, t_liste)

def normalized_entropy_over_time_animation(filename):
    # (S-S_inf)/R_fluid = 1/(gamma-1)*log(p_over_p_inf/rho_over_rho_inf**gamma)
    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename)

    gamma = dico_data["gamma"]
    M_inf = dico_data["M_inf"]
    R_fluid = dico_data["R_fluid"]
    p_0_inf = dico_data["p_0_inf"]
    T_0_inf = dico_data["T_0_inf"]

    T_inf =  T_0_inf/(1+(gamma-1)*M_inf*M_inf/2)
    p_inf = p_0_inf*(T_inf/T_0_inf)**(gamma/(gamma-1))
    rho_inf = p_inf/(R_fluid*T_inf)

    rho_all = s_all[0]
    u_all = s_all[1]/rho_all
    v_all = s_all[2]/rho_all
    U0_all = s_all[3]/rho_all
    p_all = (gamma-1)*rho_all*(U0_all - 0.5*(u_all*u_all + v_all*v_all))

    p_over_p_inf = p_all/p_inf
    rho_over_rho_inf = rho_all/rho_inf
    arg = p_over_p_inf/rho_over_rho_inf**gamma
    S_over_R_fluid = 1/(gamma-1)*np.log(arg)    


    #-----
    title_plot = "(S-S_inf)/R"
    data_name =  "Entropy "+  filename[:-4].replace("_", " ")
    data_over_time_animation(S_over_R_fluid, data_name, title_plot, x, y_all, M_inf_liste, t_liste)

def normalized_vorticity_over_time_animation(filename):
    # w*c/U_inf with w = dv/dx - du/dy
    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename)

    gamma = dico_data["gamma"]
    M_inf = dico_data["M_inf"]
    c = dico_data["c"]
    R_fluid = dico_data["R_fluid"]
    p_0_inf = dico_data["p_0_inf"]
    T_0_inf = dico_data["T_0_inf"]

    T_inf =  T_0_inf/(1+(gamma-1)*M_inf*M_inf/2)
    p_inf = p_0_inf*(T_inf/T_0_inf)**(gamma/(gamma-1))
    rho_inf = p_inf/(R_fluid*T_inf)

    rho_all = s_all[0]
    u_all = s_all[1]/rho_all
    v_all = s_all[2]/rho_all

    U_inf = p_inf / ((gamma - 1.) * rho_inf) 

    du_dy_all = np.gradient(u_all, axis=0)
    dv_dx_all = np.gradient(v_all, axis=1)

    w_all = dv_dx_all - du_dy_all

    w_normalized_all = w_all*c/U_inf

    #-----
    title_plot = "wc/U_inf"
    data_name = "Vorticity "+  filename[:-4].replace("_", " ")
    data_over_time_animation(w_normalized_all, data_name, title_plot, x, y_all, M_inf_liste, t_liste)
    