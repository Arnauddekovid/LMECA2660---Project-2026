import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os 
from scipy.ndimage import gaussian_filter  
import struct


#things t

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

def vorticity(u, v, x, y):
    """
        pre:
            u: 2D array of u speed (shape (2N, M))
            v: 2D array of v speed (shape (2N, M))
            x: 2D array of x coordinates (shape (2N, M))
            y: 2D array of y coordinates (shape( 2N, M)) 
                y[0][...] =  0+ and airfoil up
                y[N-1][...] =  upper wall
                y[N][...] =  andlower wall
                y[2N -1] = 0- and airfoil down

        post:
            returns the vorticity
    """
    N_total, M = u.shape
    N = N_total // 2  # taille d'un demi-domaine

    omega = np.zeros_like(u)

    for j in range(1, N_total - 1):

        for i in range(1, M - 1):
            du_di = 0.5 * (u[j, i+1] - u[j, i-1])
            du_dj = 0.5 * (u[j+1, i] - u[j-1, i])
            dv_di = 0.5 * (v[j, i+1] - v[j, i-1])
            dv_dj = 0.5 * (v[j+1, i] - v[j-1, i])

            dx_di = 0.5 * (x[j, i+1] - x[j, i-1])
            dx_dj = 0.5 * (x[j+1, i] - x[j-1, i])
            dy_di = 0.5 * (y[j, i+1] - y[j, i-1])
            dy_dj = 0.5 * (y[j+1, i] - y[j-1, i])

            J = dx_di * dy_dj - dx_dj * dy_di

            if abs(J) < 1e-12:
                J = 1e-12

            # dans le domaine down, dy_dj < 0 donc J < 0
            # on force J à avoir le bon signe selon le domaine
            if j >= N and J > 0:
                J = -J
            elif j < N and J < 0:
                J = -J

            du_dy = (-dy_di * du_dj + dy_dj * du_di) / J
            dv_dx = ( dx_di * dv_dj - dx_dj * dv_di) / J

            omega[j, i] = dv_dx - du_dy

    omega[0, :]  = omega[1, :]
    omega[-1, :] = omega[-2, :]
    omega[:, 0]  = omega[:, 1]
    omega[:, -1] = omega[:, -2]



    return omega

def Cd_and_Cl(dico_data, x, y, Cp):
    """
        pre:
            dico_data: dictionary containing the simulation parameters (M_inf, T_0_inf, p_0_inf, R_fluid, gamma, c, mu, omega, N, M, i_trail, ilead)
            x: 2D array of x coordinates (shape (2N, M))
            y: 2D array of y coordinates (shape (2N, M))
            Cp: 2D array of Pressure coeffcient (shape (2N, M))
        post:
            cd : the drag coefficent at the airfoil
            cl : the lift coefficent at the airfoil
    """


    # returns the drag and lift coefficient at the airfoil at time t
    i_lead = dico_data["i_lead"]
    i_trail = dico_data["i_trail"]
    c = dico_data["c"]
    
    x_airfoil = x[0][i_lead:i_trail+1]
    y_upper_airfoil = y[0][i_lead:i_trail+1]
    y_lower_airfoil = y[-1][i_lead:i_trail+1]

    Cp_upper_airfoil = Cp[0][i_lead:i_trail+1]
    Cp_lower_airfoil = Cp[-1][i_lead:i_trail+1]

    dx        = np.diff(x_airfoil)
    dy_upper  = np.diff(y_upper_airfoil)
    dy_lower  = np.diff(y_lower_airfoil)

    # upperpart
    coef =1
    dy = dy_upper
    Cp = np.diff(Cp_upper_airfoil)

    theta = np.arctan2(dy, dx) # angle of the flow at each points of the edges of the plates relative to the x axis (horizontal)
    ds = np.sqrt(dx**2 + dy**2)

    cd_up =  np.sum(coef*Cp*np.sin(theta)*ds) # the sum of the cp*ds*cos(theta) for each points of the edges of the plates
    cl_up =  np.sum(-coef*Cp*np.cos(theta)*ds) # the sum of the cp*ds*sin(theta) for each points of the edges of the plates

    #lower
    coef = -1
    dy = dy_lower
    Cp = np.diff(Cp_lower_airfoil)
    theta = np.arctan2(dy, dx) # angle of the flow at each points of the edges of the plates relative to the x axis (horizontal)
    ds = np.sqrt(dx**2 + dy**2)

    cd_down =  np.sum(coef*Cp*np.sin(theta)*ds) # the sum of the cp*ds*cos(theta) for each points of the edges of the plates
    cl_down =  np.sum(-coef*Cp*np.cos(theta)*ds) # the sum of the cp*ds*sin(theta) for each points of the edges of the plates

    cd =  -(cd_up + cd_down)
    cl =  cl_up + cl_down
    return cd, cl

def leading_and_trailing_edge_finder(x):
    """
        pre
            x: 1D array of x coordinates (shape ( M))
        post
            returns the indexes of the leading and trailing edges
    """
    i = 0
    while (x[i] < 0.):
        i +=1
    
    i_lead = i
    while (x[i] <1.):
        i +=1

    i_trail = i
    return i_lead, i_trail

def what_we_want_to_plot_at_time_t(dico_data, M_inf, x, y, s):
    """
    pre:
        dico_data: dictionary containing the simulation parameters (M_inf, T_0_inf, p_0_inf, R_fluid, gamma, c, mu, omega, N, M)
        M_inf:the Mach number at infinity
        x: 2D array of x coordinates (shape (2N, M))
        y: 3D array of y coordinates for all frames (shape(frames) 2N, M))
        s_all: 2D array of conservative variables for all frames (shape (4, 2N, M))
    post:
        Mach, Cp, total_enthalpy, total_pressure, entropy, omega
    """

    c = dico_data["c"]
    gamma = dico_data["gamma"]  
    R_fluid = dico_data["R_fluid"]
    p_0_inf = dico_data["p_0_inf"]
    T_0_inf = dico_data["T_0_inf"]

    T_inf = T_0_inf/(1+0.5*(gamma-1)*M_inf*M_inf)
    p_inf = p_0_inf*(T_inf/T_0_inf)**(gamma/(gamma-1))
    rho_inf = p_inf/R_fluid*T_inf

    u_inf = M_inf*np.sqrt(gamma*R_fluid*T_inf)
    v_inf = 0
    U_inf = T_inf*R_fluid/(gamma-1) 
    U0_inf = U_inf+u_inf*u_inf/2

    rho = s[0]
    u = s[1]/rho
    v = s[2]/rho
    U0 = s[3]/rho


    # truce en inf

    T_inf =  T_0_inf/(1+(gamma-1)*M_inf*M_inf/2)
    p_inf = p_0_inf*(T_inf/T_0_inf)**(gamma/(gamma-1))
    rho_inf = p_inf/(R_fluid*T_inf)
    u_inf = M_inf*np.sqrt(gamma*R_fluid*T_inf)

    H_0_inf = gamma/(gamma-1) * p_inf/rho_inf+0.5*u_inf*u_inf

    # trucs utile pr tout
    u_square_v_square = u *u  + v *v 
    p_over_rho  = (gamma-1)*(U0  -0.5*(u_square_v_square))

    # Mach
    Mach = np.sqrt( u_square_v_square/(gamma*p_over_rho))

    #Cp
    p = p_over_rho*rho
    num = p/p_inf-1
    den = gamma*M_inf*M_inf/2
    Cp = num/den



    # total enthalpy
    H_0 =  gamma/(gamma-1)*p_over_rho+ 0.5*u_square_v_square
    total_enthalpy = H_0/H_0_inf


    # total pressure
    p_0 = p*(1+(gamma-1)*Mach*Mach/2)**(gamma/(gamma-1))
    total_pressure = p_0/p_0_inf

    #total entropy
    p_over_p_inf = p/p_inf
    rho_over_rho_inf = rho/rho_inf
    arg = p_over_p_inf/rho_over_rho_inf**gamma
    entropy = 1/(gamma-1)*np.log(arg)    



    #vorticity
    omega = c*vorticity(u,v, x, y)/u_inf


    ### other things to add

    return [Mach, Cp, total_enthalpy, total_pressure, entropy, omega]

def get_last_frame_data_seek(filename):
    """
    pre:
        filename: path to the binary file containing the simulation data
    post:
        dico_data: some data of the grid
        x: 2D array of x coordinates (shape (2N, M))
        y_last: 2D array of the last ycoordinates (shape (2N, M))
        s_last: 2D array of s coordinates (shape (4, 2N, M))
        t_val: the last timestep
        M_inf_val: the last mach inf (not thee same as the first one if the Mach is rising
    """
    with open(filename, 'rb') as fp:
        # --- read header (fixed size) ---
        N, M = struct.unpack('ii', fp.read(8))
        keys = ['M_inf','T_0_inf','p_0_inf','R_fluid','gamma','c','mu','omega']
        vals = struct.unpack('8d', fp.read(64))
        dico_data = dict(zip(keys, vals))
        dico_data['N'] = N
        dico_data['M'] = M

        # read 1D x grid
        x1d = np.frombuffer(fp.read(M * 8), dtype=np.float64).copy()
        # expand to 2D for plotting (2N rows, M columns)
        x = np.repeat(x1d, 2*N).reshape(M, 2*N).T
        x_half = np.repeat(x1d, N).reshape(M,N).T
        # --- compute sizes ---
        header_size = 8 + 64 + M * 8
        step_size = 16 + 2 * N * M * 8 + 2 * 4 * N * M * 8
        step_size = 16 + 80 * N * M   # 80*N*M bytes per time step

        # get total file size
        fp.seek(0, os.SEEK_END)
        file_size = fp.tell()
        n_steps = (file_size - header_size) // step_size

        if n_steps <= 0:
            raise ValueError("No complete time steps found in file")

        # seek to the last step
        last_step_offset = header_size + (n_steps - 1) * step_size
        fp.seek(last_step_offset)

        # read last time step
        t_val, M_inf_val = struct.unpack('2d', fp.read(16))
        y_up    = np.frombuffer(fp.read(N*M*8), dtype=np.float64).reshape(N, M).copy()
        y_down  = np.frombuffer(fp.read(N*M*8), dtype=np.float64).reshape(N, M).copy()
        s_up    = np.frombuffer(fp.read(4*N*M*8), dtype=np.float64).reshape(4, N, M).copy()
        s_down  = np.frombuffer(fp.read(4*N*M*8), dtype=np.float64).reshape(4, N, M).copy()

        # combine up/down
        y_last = np.vstack([y_up, y_down])          # shape (2*N, M)
        s_last = np.concatenate([s_up, s_down], axis=1)  # shape (4, 2*N, M)

        return dico_data, x, y_last, s_last, t_val, M_inf_val

def all_plots_animations(filename):
    """
    pre:
        filename: path to the binary file containing the simulation data
    post:
        plots Mach, Cp, enthalpy, p0_ratio, entropy, vorticit evolution over time and save them as videos
    """

    def write_data_walls( x, datas_up,datas_down, ax, title):
        ax.clear()
        ax.plot(x, datas_up, "red", label="Upper wall")
        ax.plot(x, datas_down, "blue", label="Lower wall")
        ax.grid()
        ax.legend()
        ax.set_xlim(x[0], x[-1])
        ax.set_title(title)

    def write_data_2D(N, x, y, datas, fig, ax, vmin, vmax, color, title):

        ax.clear()

        ax.pcolormesh(
            x[N:], y[N:], datas[N:],
            cmap=color,
            vmin=vmin,
            vmax=vmax,
            shading = "gouraud",
        )

        ax.pcolormesh(
            x[N:], y[:N], datas[:N],
            cmap=color,
            vmin=vmin,
            vmax=vmax,
            shading = "gouraud",
        )

        ax.set_aspect('equal', adjustable='box')
        ax.set_ylim(-2, 2)
        if (M_inf_end >1):
            ax.set_ylim(-3, 3)
        ax.set_title(title)


    dico, x, y, s, t_end, M_inf_end = get_last_frame_data_seek(filename)

    i_lead, i_trail = leading_and_trailing_edge_finder(x[0])


    N = dico["N"]
    M = dico["M"]

    dico["i_lead"] = i_lead
    dico["i_trail"] = i_trail

    Mach, Cp, enthalpy, p0_ratio, entropy, vorticit = what_we_want_to_plot_at_time_t(dico, M_inf_end, x, y, s)



    quantities = {
        "Mach": Mach,
        "Cp": Cp,
        "total enthalpy": enthalpy,
        "total pressure": p0_ratio,
        "entropy": entropy,
        "vorticity": vorticit
    }

    vmin = {
        key: np.nanmin(val)
        for key, val in quantities.items()
    }

    vmax = {
        key: np.nanmax(val)
        for key, val in quantities.items()
    }

    vmin["Mach"] = 0.6

    colourmaps = {
        "Mach": "jet",
        "Cp": "seismic",
        "total enthalpy": "seismic",
        "total pressure": "seismic",
        "entropy": "seismic",
        "vorticity": "seismic"
    }
    
    time_values = []
    Cd = []
    Cl = []

    figures = {}
    axes = {}
    writers = {}

    names = ["Mach","Cp","total enthalpy","total pressure","entropy","vorticity"]
    names_walls = ["Mach walls", "Mach airfoil", "Cp walls","Cp airfoil"]

    readable_time = convert_time_to_readable_format()

    for name in names:
        fig, ax = plt.subplots()
        writer = animation.FFMpegWriter(fps=30)

        fig.colorbar(plt.cm.ScalarMappable(cmap=colourmaps[name], norm=plt.Normalize(vmin=vmin[name], vmax=vmax[name])), ax=ax)

        video_name = name + " " + filename.split("/")[-1].replace("_", " ")[:-4]
        string = "Animations/" + video_name+  " at " + readable_time + ".mp4"
        writer.setup(fig, string, dpi=300)
        figures[name] = fig
        axes[name] = ax
        writers[name] = writer

    for name in names_walls:
        fig, ax = plt.subplots()
        writer = animation.FFMpegWriter(fps=30)

        video_name = name + " " + filename.split("/")[-1].replace("_", " ")[:-4]
        string = "Animations/" + video_name+  " at " + readable_time + ".mp4"
        writer.setup(fig, string, dpi=300)
        figures[name] = fig
        axes[name] = ax
        writers[name] = writer



    
    ## add the other plots:

    with open(filename, 'rb') as fp:

        fp.read(8)      # N, M
        fp.read(64)     # parameters
        fp.read(M * 8)  # x

        frame_idx = 0

        while True:

            chunk = fp.read(16)

            if len(chunk) < 16:
                break

            t_val, M_inf = struct.unpack('2d', chunk)
            time_values.append(t_val)

            y_up = np.fromfile(fp,dtype=np.float64,count=N * M).reshape(N, M)
            y_down = np.fromfile(fp,dtype=np.float64,count=N * M).reshape(N, M)
            s_up = np.fromfile(fp,dtype=np.float64,count=4 * N * M).reshape(4, N, M)
            s_down = np.fromfile(fp,dtype=np.float64, count=4 * N * M).reshape(4, N, M)

            y = np.vstack([y_up, y_down])
            s = np.concatenate([s_up, s_down], axis=1)

            Mach, Cp, enthalpy, p0_ratio, entropy, vorticit = what_we_want_to_plot_at_time_t(dico,M_inf,x,y,s)
            Cd_i, Cl_i =  Cd_and_Cl(dico, x, y, Cp)

            Cd.append(Cd_i)
            Cl.append(Cl_i)

            Mach_upper_wall = Mach[N-1]
            Mach_lower_wall = Mach[N]
            Cp_upper_wall = Cp[N-1]
            Cp_lower_wall = Cp[N]

            Mach_upper_airfoil = Mach[0][i_lead:i_trail+1]
            Mach_lower_airfoil= Mach[-1][i_lead:i_trail+1]
            Cp_upper_airfoil = Cp[0][i_lead:i_trail+1]
            Cp_lower_airfoil = Cp[-1][i_lead:i_trail+1]


            data_dict = {"Mach": Mach,"Cp": Cp,"total enthalpy": enthalpy,"total pressure": p0_ratio,"entropy": entropy,"vorticity": vorticit}
            for key in names:
                title = f"{key} with M_inf = {M_inf:.4f} at t = {t_val:.6f}"
                write_data_2D(N,x,y,data_dict[key],figures[key],axes[key],vmin[key],vmax[key],colourmaps[key],title)
                writers[key].grab_frame()

            # walls stuff
            title = f"Mach at walls with M_inf = {M_inf:.4f} at t = {t_val:.6f}"
            write_data_walls(x[0], Mach_upper_wall, Mach_lower_wall, axes["Mach walls"],title)
            axes["Mach walls"].set_ylim((M_inf*(1-.2), M_inf_end*(1+.2)))
            writers["Mach walls"].grab_frame()

            title = f"Mach at airfoil with M_inf = {M_inf:.4f} at t = {t_val:.6f}"
            write_data_walls(x[0][i_lead:i_trail+1], Mach_upper_airfoil, Mach_lower_airfoil, axes["Mach airfoil"],title)
            axes["Mach airfoil"].set_ylim((M_inf_end*(1-.7), M_inf_end*(1+.7)))
            writers["Mach airfoil"].grab_frame()



            title = f"Cp at walls with M_inf = {M_inf:.4f} at t = {t_val:.6f}"
            write_data_walls(x[0], Cp_upper_wall, Cp_lower_wall, axes["Cp walls"],title)
            axes["Cp walls"].set_ylim((-0.8, 0.8))
            writers["Cp walls"].grab_frame()


            title = f"Cp at airfoil with M_inf = {M_inf:.4f} at t = {t_val:.6f}"
            write_data_walls(x[0][i_lead:i_trail+1], Cp_upper_airfoil, Cp_lower_airfoil, axes["Cp airfoil"],title)
            axes["Cp airfoil"].set_ylim((-0.8, 0.8))
            writers["Cp airfoil"].grab_frame()

            print(t_val/t_end)

            frame_idx += 1
    for key in names:
        writers[key].finish()

    fig_cd, ax_cd = plt.subplots()
    ax_cd.plot(time_values, Cd)
    ax_cd.grid()
    ax_cd.set_title("Cd evolution over time")
    ax_cd.set_xlabel("time")
    ax_cd.set_ylabel("Cd")
    plot_name = "Cd " + filename.split("/")[-1].replace("_", " ")[:-4]
    string = "Animations/" + plot_name + " at " + readable_time + ".svg"
    fig_cd.savefig(string)
    plt.close(fig_cd)

    # ── Cl plot ──────────────────────────────────────────────────────────────────
    fig_cl, ax_cl = plt.subplots()
    ax_cl.plot(time_values, Cl)
    ax_cl.grid()
    ax_cl.set_title("Cl evolution over time")
    ax_cl.set_xlabel("time")
    ax_cl.set_ylabel("Cl")
    plot_name = "Cl " + filename.split("/")[-1].replace("_", " ")[:-4]
    string = "Animations/" + plot_name + " at " + readable_time + ".svg"
    fig_cl.savefig(string)
    plt.close(fig_cl)
    print("\nVideos successfully created.")

# plots the data by extracting all the frames and then creating an animation with matplotlib, useful to debug

def extract_data(filename):
    """
        pre:
            filename: path to the binary file containing the simulation data
        post:
            dico_data: dictionary containing the simulation parameters (M_inf, T_0_inf, p_0_inf, R_fluid, gamma, c, mu, omega, N, M)
            t_list: list of time values for each frame
            M_inf_list: list of Mach numbers for each frame
            x: 2D array of x coordinates (shape (2N, M))
            y_all: 3D array of y coordinates for all frames (shape (n_frames, 2N, M))
            s_all: 4D array of conservative variables for all frames (shape (4, n_frames, 2N, M))
    """

    dico_data = {}
    t_list    = []
    M_inf_list = []
    y_frames  = []
    s_frames  = []

    with open(filename, 'rb') as fp:

        # ── header ints ──────────────────────────────────────────────────────
        N, M = struct.unpack('ii', fp.read(8))

        # ── header doubles ───────────────────────────────────────────────────
        keys = ['M_inf','T_0_inf','p_0_inf','R_fluid','gamma','c','mu','omega']
        vals = struct.unpack('8d', fp.read(64))
        dico_data = dict(zip(keys, vals))
        dico_data['N'] = N
        dico_data['M'] = M

        # ── x grid ───────────────────────────────────────────────────────────
        x = np.frombuffer(fp.read(M * 8), dtype=np.float64).copy()
        x = np.repeat(x,2*N).reshape(M,2*N).T

        # ── time steps ───────────────────────────────────────────────────────
        frame_size = (2 + 2*N + 2*4*N) * M * 8 + 16
        # = t + M_inf (16 bytes) + y_up + y_down + s_up + s_down

        while True:
            chunk = fp.read(16)
            if len(chunk) == 0:
                break
            if len(chunk) < 16:
                print("Warning: truncated frame, stopping.")
                break

            t_val, M_inf_val = struct.unpack('2d', chunk)
            t_list.append(t_val)
            M_inf_list.append(M_inf_val)

            # y_up shape (N, M)
            y_up   = np.frombuffer(fp.read(N*M*8), dtype=np.float64).reshape(N, M).copy()
            # y_down shape (N, M)
            y_down = np.frombuffer(fp.read(N*M*8), dtype=np.float64).reshape(N, M).copy()

            # s_up shape (4, N, M)
            s_up   = np.frombuffer(fp.read(4*N*M*8), dtype=np.float64).reshape(4, N, M).copy()
            # s_down shape (4, N, M)
            s_down = np.frombuffer(fp.read(4*N*M*8), dtype=np.float64).reshape(4, N, M).copy()

            # stack vertically: up domain on top, down domain below
            y_full = np.vstack([y_up, y_down])            # (2N, M)
            
            s_full = np.concatenate([s_up, s_down], axis=1)  # (4, 2N, M)

            y_frames.append(y_full)
            s_frames.append(s_full)

    y_all = np.stack(y_frames, axis=0)      # (n_frames, 2N, M)
    s_all = np.stack(s_frames, axis=0)      # (n_frames, 4, 2N, M)
    # reorder to (4, n_frames, 2N, M) to match your existing code
    s_all = np.moveaxis(s_all, 1, 0)

    return dico_data, np.array(t_list), np.array(M_inf_list), x, y_all, s_all

def evolution_animation_grid(filename):
    """
        pre:
            filename: path to the binary file containing the simulation data
        post:
            plots the grid evolution over time by extracting all the frames
    """

    dico_data, t_liste,_, x, y_all, _= extract_data(filename)
    N = dico_data["N"]
    M = dico_data["M"]

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

    string = "Animations/grid "+filename.split("/")[-1].replace("_", " ")[:-4] + " at " + time + ".mp4"

    anim.save(string, writer='ffmpeg', fps=30)  

def data_over_time_animation(datas, video_name, plot_title,color_style,  x, y_all, M_inf_liste, t_liste):
    """
        pre:
            data : the data we want to plot over time
            video_name : the name of the video to be created (without extension) eg "vorticity low res M inf 0.85 oscillating mu 0.05 until 0.01s "
            plot_title : the title to be displayed on the plot eg "Vorticity "
            x : 2D array of x coordinates (shape (2N, M))
            y_all : 3D array of y coordinates for all frames (shape (n_frames, 2N, M))
            M_inf_liste : list of Mach numbers for each frame
            t_liste : list of time values for each frame
        post:
            plots the grid evolution over time by extracting all the frames
    """
    vmin = min(np.min(data) for data in datas)
    vmax = max(np.max(data) for data in datas)
    N = len(x)//2
    fig, ax = plt.subplots()

    M_inf = M_inf_liste[0]


    mesh_up = ax.pcolormesh(
        x[N:],
        y_all[0][N:],
        datas[0][N:],
        vmin=vmin,
        vmax=vmax,            
        shading = "gouraud", 
        cmap=color_style
    )

    mesh_down = ax.pcolormesh(
        x[N:],
        y_all[0][:N],
        datas[0][:N],
        vmin=vmin,
        vmax=vmax,
        shading = "gouraud",
        cmap=color_style

    )

    cbar = plt.colorbar(mesh_up, ax=ax)
    title = ax.set_title(f"{plot_title} with M_inf = {M_inf_liste[0]:.3f} t = {t_liste[0]:.8f}")
    ax.axis('equal')

    # --- Update function ---
    def update(frame):
        ax.clear()

        print(frame/len(t_liste ))

        mesh_up = ax.pcolormesh(
            x[N:],
            y_all[frame][N:],
            datas[frame][N:],
            vmin=vmin,
            vmax=vmax,
            shading = "gouraud",
            cmap=color_style
        )

        mesh_down = ax.pcolormesh(
            x[N:],
            y_all[frame][:N],
            datas[frame][:N],
            vmin=vmin,
            vmax=vmax,
            shading = "gouraud",
            cmap=color_style

        )

        ax.set_ylim(-2, 2)
        if M_inf >1:
            ax.set_ylim(-3, 3)

        # Puis impose le ratio égal sans toucher aux limites
        ax.set_aspect('equal', adjustable='box')

        ax.set_title(f"{plot_title}with M_inf = {M_inf_liste[frame]:.3f} t = {t_liste[frame]:.8f}")

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
    string = "Animations/" + video_name+  " at " + time + ".mp4"

    anim.save(string, writer='ffmpeg', fps=30)  

def Mach_number_over_time_animation(filename):
    """
    pre:
        filename: path to the binary file containing the simulation data
    post:
        plots the Mach number evolution over time by extracting all the frames and save it as a video
    """
    dico_data, t_liste,M_inf_liste, x, y_all, s_all= extract_data(filename,)

    gamma = dico_data["gamma"]
    rho_all = s_all[0]
    u_all = s_all[1]/rho_all
    v_all = s_all[2]/rho_all
    U0_all = s_all[3]/rho_all

    p_all_over_rho_all = (gamma-1)*(U0_all -0.5*(u_all*u_all + v_all*v_all))

    Mach_all = np.sqrt( (u_all*u_all + v_all*v_all)/(gamma*p_all_over_rho_all))


    #-----
    title_plot = "Mach number "
    video_name = title_plot+filename.split("/")[-1].replace("_", " ")[:-4]
    color_style = 'jet'
    data_over_time_animation(Mach_all, video_name, title_plot,color_style,  x, y_all, M_inf_liste, t_liste)

def Cp_over_time_animation(filename):
    """
    pre:
        filename: path to the binary file containing the simulation data
    post:
        plots the Cp evolution over time by extracting all the frames and save it as a video
    """
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
    title_plot = "Cp "
    video_name = title_plot+filename.split("/")[-1].replace("_", " ")[:-4]
    color_style = 'jet'
    data_over_time_animation(Cp_all, video_name, title_plot,color_style,  x, y_all, M_inf_liste, t_liste)

def normalized_total_enthalpy_over_time_animation(filename):
    """
    pre:
        filename: path to the binary file containing the simulation data
    post:
        plots the normalized total enthalpy evolution over time by extracting all the frames and save it as a video
    """
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
    title_plot = "H_0/H_0_inf "
    video_name = "enthalpy "+filename.split("/")[-1].replace("_", " ")[:-4]
    color_style = 'seismic'
    data_over_time_animation(H_0_over_H_0_inf, video_name, title_plot,color_style,  x, y_all, M_inf_liste, t_liste)

def normalized_total_pressure_over_time_animation(filename):
    """
    pre:
        filename: path to the binary file containing the simulation data
    post:
        plots the normalized total pressure evolution over time by extracting all the frames and save it as a video
    """
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

    title_plot = "p_0/p_0_inf "
    video_name =   "total pressure "+filename.split("/")[-1].replace("_", " ")[:-4]
    color_style = 'seismic'
    data_over_time_animation(p_O_over_p_0_inf, video_name, title_plot,color_style,  x, y_all, M_inf_liste, t_liste)

def normalized_entropy_over_time_animation(filename):
    """
    pre:
        filename: path to the binary file containing the simulation data
    post:
        plots the normalized total entropy evolution over time by extracting all the frames and save it as a video
    """
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
    data_name =  "Entropy "+filename.split("/")[-1].replace("_", " ")[:-4]
    color_style = 'seismic'
    data_over_time_animation(S_over_R_fluid, data_name, title_plot,color_style,  x, y_all, M_inf_liste, t_liste)

def normalized_vorticity_over_time_animation(filename):
    """
    pre:
        filename: path to the binary file containing the simulation data
    post:
        plots the normalized vorticity evolution over time by extracting all the frames and save it as a video
    """
    dico_data, t_liste, M_inf_liste, x, y_all, s_all = extract_data(filename)

    gamma = dico_data["gamma"]
    M_inf = dico_data["M_inf"]
    c = dico_data["c"]
    R_fluid = dico_data["R_fluid"]
    p_0_inf = dico_data["p_0_inf"]
    T_0_inf = dico_data["T_0_inf"]

    T_inf = T_0_inf / (1.0 + (gamma - 1.0) * M_inf * M_inf / 2.0)
    p_inf = p_0_inf * (T_inf / T_0_inf) ** (gamma / (gamma - 1.0))
    rho_inf = p_inf / (R_fluid * T_inf)

    rho_all = s_all[0]
    u_all = s_all[1] / rho_all
    v_all = s_all[2] / rho_all

    # correct free-stream velocity
    U_inf = M_inf*np.sqrt(R_fluid*gamma*T_inf)

    vorticity_array = []

    for t in range(len(t_liste)):
        vorticity_t = vorticity(u_all[t], v_all[t], x, y_all[t])
        vorticity_array.append(vorticity_t)

    vorticity_all = np.array(vorticity_array)*c/U_inf
    #-----
    title_plot = "wc/U_inf"
    data_name = "Vorticity "+ filename.split("/")[-1].replace("_", " ")[:-4]
    color_style = 'seismic'
    data_over_time_animation(vorticity_all, data_name, title_plot,color_style, x, y_all, M_inf_liste, t_liste)
