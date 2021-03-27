
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

def Ez_analytical(x,y,a=1.0,k=2*np.pi, rho_d=3.0, phi_d=np.pi, S=20):

    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)

    E0 = 0
    E1 = 0

    for n in range(-S,S+1):

        c = - jv(n, k*a)/hankel2(n, k*a)

        E0 += hankel2(n, k*rho_d)*(jv(n, k*rho) + c*hankel2(n, k*rho))*np.exp(1j*n*(phi-phi_d))
        E1 += hankel2(n, k*rho)*(jv(n, k*rho_d) + c*hankel2(n, k*rho_d))*np.exp(1j*n*(phi-phi_d))

    mask = rho <= rho_d

    E1[mask] = E0[mask]

    return np.real(E1)

def get_scale(alpha, a, b):
    return np.sqrt(np.sum((alpha*a-b)**2))

def check():

    plt.rc("font", family="serif")
    plt.rc("font", size=8)
    plt.rc("mathtext", fontset="cm")
    # matplotlib.rc('text', usetex = True)
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    

    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[0.35, 1], wspace=0.05, bottom=0.14, top=0.97, left=0.1, right=0.85)

    

    #==============================================================================
    # Simulation results
    #==============================================================================

    # get a list of all the files in this directory
    files = get_files('.', include=['plt'], exclude=["temp"], get_all=True)

    f = files[-1]

    fig = plt.figure(figsize=(5,3))
    ax1 = fig.add_subplot(gs[0,1])

    data = ReadBoxLib(f)
    time = data.time
        
    data = ReadBoxLib(f, max_level=-1)

    xn, Dz = data.get("z_D-field", grid="node")
    xc, vf = data.get("vfrac-field")
    data.close()

    yn, xn = np.meshgrid(xn[1], xn[0])
    yc, xc = np.meshgrid(xc[1], xc[0])

    # maske out the cylinder
    r = xc**2 + yc**2
    mask = r <= 1
    xc = np.ma.masked_where(mask,xc)
    yc = np.ma.masked_where(mask,yc)

    # plot field
    nx, ny = xn.shape
    ii = int(nx/2)

    merge = np.zeros(xc.shape)

    # upper half

    slicer = (slice(0,nx), slice(int(ny/2),ny))
    Dz = Dz[slicer]
    merge[slicer] = Dz

    # analytical solution
    slicer = (slice(0,nx), slice(0,int(ny/2)))
    

    try:
        Ez = pickle.load(open("Ez.p", "rb"))
    except:
        Ez = Ez_analytical(xc[slicer],yc[slicer],a=1.0,k=2*np.pi, rho_d=3.0, phi_d=np.pi, S=20)
        pickle.dump(Ez, open("Ez.p", "wb"))

    # scale
    D_sample = Dz[ii,:]
    E_sample = Ez[ii,:]

    alpha = D_sample.max()/E_sample.max()

    print("scaled by alpha = ",alpha)

    # alpha = fmin(get_scale, 1.0, (E_sample, D_sample))

    merge[slicer] = -alpha*Ez
    E_sample *= -alpha

    #L2 = np.sqrt(np.sum(E_sample - D_sample))

    # blank out geometry

    cm = ax1.pcolormesh(xn, yn, merge)
    cs = ax1.contour(xc, yc, vf, levels=[0.5], linewidths=[0.5], colors=['k'])

    for line in cs.allsegs[0]:

        x = line[:,0].tolist()
        y = line[:,1].tolist()

        plt.fill(x, y, 'w', hatch='////')

    
    ax1.set_aspect(1)
    ax1.set_yticklabels([])
    ax1.set_xlabel(r"$x$")
    ax1_divider = make_axes_locatable(ax1)
    # add an axes to the right of the main axes.
    cax = ax1_divider.append_axes("right", size="7%", pad="2%")
    cb = plt.colorbar(cm, cax=cax, label=r"$D_z$")

    ax1.text(0.95, 0.95, r"sim.", horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes)

    ax1.text(0.95, 0.05, r"ref.", horizontalalignment='right', verticalalignment='bottom', transform=ax1.transAxes)



    ax2 = fig.add_subplot(gs[0,0])

    ii = int(nx/2)
    x = xc[ii,0]
    y = yc[ii,:]

    ax1.plot(y.size*[x], y, "r--", lw=1)

    x = xc[:,0]
    x = np.ma.masked_where(np.abs(x) < 1, x)
    ax1.plot(x, x.size*[0], 'k--', lw=0.5)

    ax2.plot(D_sample, y[0:int(ny/2)][::-1], 'k', lw=0.5)
    ax2.plot(E_sample, y[0:int(ny/2)], 'k--', lw=0.5)

    rms = np.sqrt(np.sum((E_sample - D_sample[::-1])**2)/float(len(E_sample)))

    #print("error = ",rms)


    ax2.plot(D_sample, y[int(ny/2)::], 'k', lw=0.5)
    ax2.plot(E_sample, y[int(ny/2)::][::-1], 'k--', lw=0.5)

    lim = ax2.get_xlim()
    ax2.fill_between(lim, [1,1], [-1,-1], facecolor='none', edgecolor="k", hatch='////', linewidth=0.25)

    ax2.set_xlim(lim)
    ax2.set_ylim(-8,8)
    ax2.set_xlabel(r"$D_z$")
    ax2.set_ylabel(r"$y$")

    # fig.tight_layout()

    fig.savefig("plot.png",dpi=300)
    plt.close(fig)

    if rms > 1e-3:
        return 1
    else:
        return 0


if __name__ == "__main__":
    sys.exit(check())
