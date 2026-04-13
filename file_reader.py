import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

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
        x_str = file[0].split(" ")[:-1]
        t = float(file[1])

        for x_j in x_str:
            x.append(float(x_j))
        i = 2
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

            y_all.append(y)
            u_all.append(u)

    



        """


        #plt.scatter(x,y, s=1)
        #plt.pcolormesh(x[begin:end], y[begin:end],u, edgecolors='black', linewidth=0.2)
        #plt.pcolormesh(x[begin:end], y[begin:end],np.zeros_like(x), facecolor='none', edgecolors='black', linewidth=0.5)

        plt.pcolormesh(x, y,np.zeros_like(x), facecolor='none', edgecolors='black', linewidth=0.5)
        #plt.pcolormesh(x, y, u, facecolor='none', linewidth=0.5)

        #plt.grid()
        plt.colorbar()
        plt.show()
        """
    np.array(t_liste)
    np.array(y_all)
    np.array(u_all)

    return t_liste, x, y_all, u_all

def plot_D_data(filename):
    x = []
    with open(filename, 'r') as f:

        file = f.readlines()
        x_str = file[0].split(" ")[:-1]
        t = float(file[1])
        for x_j in x_str:
            x.append(float(x_j))
        i = 2
        np.array(x)
        M = 0
        N = 0

        y = []
        D = []

        while (len(file[i].split(" "))!=2):
            y_i = file[i].split(" ")
            y_liste = []
            for y_ij in y_i[:-1]:
                y_liste.append(float(y_ij))
            y.append(np.array(y_liste))
            i+=1
        
        i+=1

        while (len(file[i].split(" "))!=2):
            D_i = file[i].split(" ")
            D_liste = []
            for D_ij in D_i[:-1]:

                D_liste.append(float(D_ij))
            D.append(np.array(D_liste))
            i+=1 

    np.array(y_liste)

    M = len(x)
    N = len(y)
    x = np.repeat(x, len(y)).reshape(M,N).T
    
    y = np.reshape(y, (N, M))

    np.array(D)    

    #plt.scatter(x,y, s=1)
    #plt.pcolormesh(x[begin:end], y[begin:end],u, edgecolors='black', linewidth=0.2)
    #plt.pcolormesh(x[begin:end], y[begin:end],np.zeros_like(x), facecolor='none', edgecolors='black', linewidth=0.5)

    plt.pcolormesh(x, y,D, edgecolors='black', linewidth=0.5)
    #plt.pcolormesh(x, y, u, facecolor='none', linewidth=0.5)

    #plt.grid()
    plt.colorbar()
    plt.show()
    
def evolution_animation(filename, show):
    """
    pre :
        filename the file we must plot 
    post:
        plots the animation of the temporal evoloution of the data described in the file filename. the data is extracted from the file using extract data
    """

    t, x, y_all, u_all = extract_data(filename)

    fig, ax = plt.subplots()
    vmin = min(np.min(u) for u in u_all)
    vmax = max(np.max(u) for u in u_all)
    extent = [x.min(), x.max(), y_all[0].min(), y_all[0].max()]

    im = ax.imshow(u_all[0], origin='lower', extent=extent, aspect='auto', cmap='viridis', vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax)
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    def update(frame_idx):
        t_i = t[frame_idx]
        im.set_data(u_all[frame_idx])
        ax.set_title(f'M_inf = {filename[8:12]} : t = {t_i}')
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=len(t),interval=10, blit=False, repeat=True)
    time = convert_time_to_readable_format() 

    string = "Animations/"+filename[:-4] + " at " + time + ".mp4"

    ani.save(string, writer='ffmpeg', fps=30)  

    if show:
        plt.show()
