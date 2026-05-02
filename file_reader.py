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

    x = []
    t_liste = []
    y_all = []
    u_all = []


    with open(filename, 'r') as f:

        file = f.readlines()
        infos_str = file[0][:-1]
        x_str = file[1].split(" ")[:-1]
        t = float(file[2])

        for x_j in x_str:
            x.append(float(x_j))
        i = 3
        np.array(x)
        M = 0
        N = 0

        x_reshape_already_done = False


        while file[i-1].split(" ")[0] != "end" :

            y = []
            u = []

            t_i = file[i-1].split(" ")[0]
            t_liste.append(float(t_i))

            while (len(file[i].split(" "))!=2):
                y_i = file[i].split(" ")
                y_liste = []
                for y_ij in y_i[:-1]:

                    y_liste.append(float(y_ij))
                y.append(np.array(y_liste))
                i+=1
            
            i+=1

            while (len(file[i].split(" "))!=2):
                u_i = file[i].split(" ")
                u_liste = []
                for u_ij in u_i[:-1]:

                    u_liste.append(float(u_ij))
                u.append(np.array(u_liste))
                i+=1 

            i+=1

            np.array(y)
            np.array(u)

            if not x_reshape_already_done:
                M = len(x)
                N = len(y)
                x = np.repeat(x, len(y)).reshape(M,N).T
                x_reshape_already_done = True
            
            y = np.reshape(y, (N, M))
            u = np.reshape(u, (N, M))   # <-- ADD THIS LINE


            y_all.append(y)
            u_all.append(u)

    np.array(t_liste)
    np.array(y_all)
    np.array(u_all)

    return infos_str, t_liste, x, y_all, u_all

def plot_other_data(filename, data_name):
    x = []
    with open(filename, 'r') as f:

        file = f.readlines()
        infos_str = file[0][:-1]
        x_str = file[1].split(" ")[:-1]
        t = float(file[2])
        for x_j in x_str:
            x.append(float(x_j))
        i = 3
        np.array(x)
        M = 0
        N = 0

        y = []
        data_liste_all = []

        while (len(file[i].split(" "))!=2):
            y_i = file[i].split(" ")
            y_liste = []
            for y_ij in y_i[:-1]:
                y_liste.append(float(y_ij))
            y.append(np.array(y_liste))
            i+=1
        
        i+=1

        while (len(file[i].split(" "))!=2):
            data_i = file[i].split(" ")
            data_liste= []
            for data_ij in data_i[:-1]:

                data_liste.append(float(data_ij))
            data_liste_all.append(np.array(data_liste))
            i+=1 

    np.array(y_liste)

    M = len(x)
    N = len(y)
    x = np.repeat(x, len(y)).reshape(M,N).T
    
    y = np.reshape(y, (N, M))

    np.array(data_liste_all)    

    #plt.scatter(x,y, s=1)
    #plt.pcolormesh(x[begin:end], y[begin:end],u, edgecolors='black', linewidth=0.2)
    #plt.pcolormesh(x[begin:end], y[begin:end],np.zeros_like(x), facecolor='none', edgecolors='black', linewidth=0.5)



    plt.pcolormesh(x[:N//2], y[:N//2] ,data_liste_all[:N//2], edgecolors='black',shading = "gouraud", linewidth=0.5)
    plt.pcolormesh(x[:N//2], y[N//2:] ,data_liste_all[N//2:], edgecolors='black',shading = "gouraud", linewidth=0.5)

    #plt.pcolormesh(x, y, u, facecolor='none', linewidth=0.5)

    #plt.grid()
    plt.title("visu of the data from "+ data_name)
    plt.colorbar()
    plt.show()

def plot_grid(filename):
    infos_str, _, x, y, _ = extract_data(filename)
    N = len(x)//2
    x_half = x[N:]
    y_up = y[0][:N]
    y_down = y[0][N:]
    plt.pcolormesh(x_half, y_up,np.zeros_like(x_half), facecolor='none',edgecolors='black', linewidth=0.2)
    plt.pcolormesh(x_half, y_down,np.zeros_like(x_half), facecolor='none',edgecolors='black', linewidth=0.2)
    plt.show()

def evolution_animation(filename, delete_txt_file_after_animation):

    # --- Load data ---
    infos_str, t_liste, x, y_all, u_all = extract_data(filename)
    N = len(x)//2


    # --- Global color scale ---
    vmin = min(np.min(u) for u in u_all)
    vmax = max(np.max(u) for u in u_all)

    # --- Figure setup ---
    fig, ax = plt.subplots()

    mesh_up = ax.pcolormesh(
        x[N:],
        y_all[0][N:],
        u_all[0][N:],
        vmin=vmin,
        vmax=vmax,            
        shading = "gouraud"

    )

    mesh_down = ax.pcolormesh(
        x[N:],
        y_all[0][:N],
        u_all[0][:N],
        vmin=vmin,
        vmax=vmax,
        shading = "gouraud"

    )

    cbar = plt.colorbar(mesh_up, ax=ax)
    title = ax.set_title(f"t = {t_liste[0]:.3f}")

    # --- Update function ---
    def update(frame):
        ax.clear()

        mesh_up = ax.pcolormesh(
            x[N:],
            y_all[frame][N:],
            u_all[frame][N:],
            vmin=vmin,
            vmax=vmax,
            shading = "gouraud"
        )

        mesh_down = ax.pcolormesh(
            x[N:],
            y_all[frame][:N],
            u_all[frame][:N],
            vmin=vmin,
            vmax=vmax,
            shading = "gouraud"

        )

        ax.set_title(infos_str + filename[:11]+f" t = {t_liste[frame]:.6f}")

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

    string = "Animations/"+filename[:-4] + " at " + time + ".mp4"

    anim.save(string, writer='ffmpeg', fps=30)  


    if delete_txt_file_after_animation:
        os.remove(filename)
    
