
import os, sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
from scipy.special import jv, hankel2
import pylab as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import pickle
from scipy.optimize import fmin

def check():

    plt.rc("font", family="serif")
    plt.rc("font", size=8)
    plt.rc("mathtext", fontset="cm")
    # matplotlib.rc('text', usetex = True)
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    

    gs = gridspec.GridSpec(nrows=3, ncols=3) #, width_ratios=[0.35, 1], wspace=0.05, bottom=0.14, top=0.97, left=0.1, right=0.85)

    

    #==============================================================================
    # Simulation results
    #==============================================================================

    # get a list of all the files in this directory
    files = get_files('.', include=['plt'], exclude=["temp"], get_all=True)

    f = files[-1]

    fig = plt.figure(figsize=(8,5))

    data = ReadBoxLib(f)
    time = data.time
        
    data = ReadBoxLib(f, max_level=-1, limits=[[-5,-5],[8,5]])

    # ---
    xn, rho_e = data.get("rho-electron1", grid="node")
    xn, rho_i = data.get("rho-ion1", grid="node")

    xn, m_e = data.get("mass-electron1", grid="node")
    xn, m_i = data.get("mass-ion1", grid="node")

    xn, q_e = data.get("charge-electron1", grid="node")
    xn, q_i = data.get("charge-ion1", grid="node")

    nd_e1 = rho_e/m_e
    nd_i1 = rho_i/m_i

    cd1 = (q_e*rho_e/m_e + q_i*rho_i/m_i)

    # ---

    xn, rho_e = data.get("rho-electron2", grid="node")
    xn, rho_i = data.get("rho-ion2", grid="node")

    xn, m_e = data.get("mass-electron2", grid="node")
    xn, m_i = data.get("mass-ion2", grid="node")

    xn, q_e = data.get("charge-electron2", grid="node")
    xn, q_i = data.get("charge-ion2", grid="node")

    nd_e2 = rho_e/m_e
    nd_i2 = rho_i/m_i

    cd2 = (q_e*rho_e/m_e + q_i*rho_i/m_i)

    # ---

    xn, Bx = data.get("x_B-field1", grid="node")
    xn, By = data.get("y_B-field1", grid="node")
    xn, Bz = data.get("z_B-field1", grid="node")

    B1 = np.sqrt(Bx**2 + By**2 + Bz**2)

    xn, Dx = data.get("x_D-field1", grid="node")
    xn, Dy = data.get("y_D-field1", grid="node")
    xn, Dz = data.get("z_D-field1", grid="node")

    D1 = np.sqrt(Dx**2 + Dy**2 + Dz**2)

    # ---

    xn, Bx = data.get("x_B-field2", grid="node")
    xn, By = data.get("y_B-field2", grid="node")
    xn, Bz = data.get("z_B-field2", grid="node")

    B2 = np.sqrt(Bx**2 + By**2 + Bz**2)

    xn, Dx = data.get("x_D-field2", grid="node")
    xn, Dy = data.get("y_D-field2", grid="node")
    xn, Dz = data.get("z_D-field2", grid="node")

    D2 = np.sqrt(Dx**2 + Dy**2 + Dz**2)

    # ---

    xc, vf_fluid = data.get("vfrac-ion1")
    xc, vf_field = data.get("vfrac-field1")
    data.close()

    yn, xn = np.meshgrid(xn[1], xn[0])
    yc, xc = np.meshgrid(xc[1], xc[0])

    plot = [
        {"data":[cd1,cd2], "label":r"$\varrho_c$", "loc":[0,0], "eb":vf_fluid, "cmap":"bwr"},
        {"data":[B1, B2], "label":r"$\left|\mathbf{B}\right|$", "loc":[1,0], "eb":vf_field, "cmap":"viridis"},
        {"data":[D1,D2], "label":r"$\left|\mathbf{D}\right|$", "loc":[2,0], "eb":vf_field, "cmap":"viridis"},
        

        {"data":[nd_e1,nd_e2], "label":r"$n_e$", "loc":[0,1], "eb":vf_fluid, "cmap":"viridis"},
        {"data":[nd_i1,nd_i2], "label":r"$n_i$", "loc":[1,1], "eb":vf_fluid, "cmap":"viridis"},
        ]

    axes = []
    
    for p in plot:

        ax = fig.add_subplot(gs[p["loc"][0], p["loc"][1]]); axes.append(ax)

        ni, nj = xc.shape

        plot_data = np.hstack((p["data"][0][:,0:int(nj/2)], p["data"][1][:,int(nj/2)::]))

        color = "w"
        if p["cmap"] == "bwr":
            big = np.max(np.abs(plot_data))
            vmin = -big
            vmax =  big
            color="k"
        else:
            vmin = np.min(plot_data)
            vmax = np.max(plot_data)

        cm = ax.pcolormesh(xn, yn, plot_data, cmap=p["cmap"], vmin=vmin, vmax=vmax)



        cs = ax.contour(xc, yc, p["eb"], levels=[0.5], linewidths=[0.5], colors=['k'])

        for line in cs.allsegs[0]:

            x = line[:,0].tolist()
            y = line[:,1].tolist()

            ax.fill(x, y, 'w', hatch='////////')

        ax.plot(ax.get_xlim(), 2*[0], "--", lw=0.5, color=color)

        
        ax.set_aspect(1)
        ax.set_yticks([])
        ax.set_xticks([])
        ax_divider = make_axes_locatable(ax)
        # add an axes to the right of the main axes.
        cax = ax_divider.append_axes("right", size="7%", pad="2%")
        cb = plt.colorbar(cm, cax=cax, label=p["label"])

    fig.savefig("plot.png",dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(check())