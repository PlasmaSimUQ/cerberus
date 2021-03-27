
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
import matplotlib.ticker as ticker

def check():

    #==============================================================================
    # Simulation results
    #==============================================================================

    # get a list of all the files in this directory
    files = get_files('.', include=['plt'], exclude=["temp"], get_all=True)

    f = files[-1]

    data = ReadBoxLib(f)
    t = data.time
        
    data = ReadBoxLib(f, max_level=-1)

    xc, u = data.get("x_vel-air")
    xc, v = data.get("y_vel-air")

    vel = np.sqrt(u**2 + v**2)

    yc, xc = np.meshgrid(xc[1], xc[0])

    R = np.sqrt(xc**2 + yc**2)

    R_linear = np.ravel(R)
    vel_linear = np.ravel(vel)

    r_max = 8.0
    R_linear = np.ma.masked_where(R_linear>r_max, R_linear)
    vel_linear = np.ma.masked_where(R_linear>r_max, vel_linear)

    I = np.argsort(R_linear)
    R_linear = R_linear[I]
    vel_linear = vel_linear[I]

    # =============================================================================
    # analytical solution
    # =============================================================================

    # D. J. Munoz, V. Springel, R. Marcus, M. Vogelsberger, L. Hernquist, 
    # Multidimensional, compressible viscous flow on a moving Voronoi mesh, 
    # Monthly Notices of the Royal Astronomical Society, 
    # Volume 428, Issue 1, 1 January 2013, Pages 254-279, 
    # https://doi.org/10.1093/mnras/sts015

    G = 1.0
    mu0 = 0.08
    rho0 = 1.0
    nu = mu0/rho0
    t0 = 10.0
    def vtheta(R,t):
        return G/(2*np.pi*R)*(1-np.exp(-R**2/(4*nu*t)))

    vt = vtheta(R_linear, data.time+t0)

    # =============================================================================
    # check
    # =============================================================================

    success = 0

    rel_err = np.abs((vel_linear - vt)/vt)

    if np.max(rel_err) > 0.01:
        success = 1
    
    # =============================================================================
    # plot
    # =============================================================================

    plt.rc("font", family="serif")
    plt.rc("font", size=8)
    plt.rc("mathtext", fontset="cm")
    # matplotlib.rc('text', usetex = True)
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    fig = plt.figure(figsize=(5,2))

    ax = fig.add_subplot(111)
    ax.plot(R_linear, vel_linear,'.', ms=2, mfc='none')
    ax.plot(R_linear, vt, 'k--', lw=1)
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$v_\theta$")

    ax = ax.twinx()

    ax.plot(R_linear, rel_err*1000, 'r.', ms=0.5)
    ax.set_ylabel(r'$\left| \frac{\hat{v}_\theta - v_\theta}{v_\theta} \right|\times 10^3$')

    ax.set_xlim(0,8)

    ylim = ax.get_ylim()
    ax.set_ylim(0, ylim[1])

    fig.tight_layout()
    fig.savefig("plot.pdf", dpi=300)

    return success
        

if __name__ == "__main__":
    sys.exit(check())