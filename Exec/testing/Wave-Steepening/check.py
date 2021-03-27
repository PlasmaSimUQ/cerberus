
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files
from LandauSimpleWave import LandauWave

import numpy as np
import pylab as plt
import matplotlib.ticker as ticker

import subprocess, shlex

def get_name(levels):
    return "Landau_n%i.plt"%levels

def run(d):

    inputs = "Landau.inputs"
    EXE = "../../local/MFP.1d.gnu.MPI.ex"

    ncells = d["n"]
    name = get_name(ncells)

    d["name"] = name
    t = d["t"]

    if 1:
        run = "rm -r {base_name}*; mpirun -n 1 {exe} {input} amr.n_cell={n} amr.plot_file={base_name} stop_time={time} 2>&1 | tee {base_name}.log".format(n=ncells, exe=EXE, base_name=name, time=t, input=inputs)
        subprocess.call(run, shell=True)

    return

def RMS(a, b):
    return np.sqrt(np.sum((a-b)**2))/float(len(a))

def get_data(f):

    rh5 = ReadBoxLib(f)
    t = rh5.time
        
    rh5 = ReadBoxLib(f, max_level=-1)

    name = "neutral"

    xc, rho = rh5.get("rho-%s"%name)
    xn, u = rh5.get("x_vel-%s"%name, grid='node')
    _, T = rh5.get("T-%s"%name)

    rh5.close()

    return t, xc[0], rho, u, T


def check():

    u0 = 0.1
    p0 = 1.0
    rho0 = 1.0
    gam = 5.0/3.0
    T0 = p0/rho0

    tbreak = np.abs(2/(u0*np.pi*(gam + 1)))

    data = [
        {"n":16},
        {"n":32},
        {"n":64},
        {"n":128},
    ]

    for d in data:

        d["t"] = 0.8*tbreak

        run(d)


        # get a list of all the files in this directory
        files = get_files('.', include=[d["name"], '.plt'], exclude=["temp", "log"], get_all=True)

        t, x, rho, u, T = get_data(files[-1])

        d["sim"] = [x, rho, u, T]
            

        # =============================================================================
        # analytical solution
        # =============================================================================

        LW = LandauWave(rho0, T0, u0, gam, x=x)
        LW.getSol(t)

        d["ref"] = [LW.x, LW.rho, LW.u, LW.T]

        # ERROR
        d["RMS"] = [RMS(LW.rho, rho), RMS(LW.u, u), RMS(LW.T, T)]


    # =============================================================================
    # plot
    # =============================================================================

    plt.rc("font", family="serif")
    plt.rc("font", size=8)
    plt.rc("mathtext", fontset="cm")
    # matplotlib.rc('text', usetex = True)
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    plt.rcParams.update(params)


    axes = []
    fig = plt.figure(figsize=(5,3.5))

    nr = 3
    nc = 2

    x, rho, u, T = data[-1]["sim"]
    x, r_rho, r_u, r_T = data[-1]["ref"]

    N = []
    err = []

    for d in data:
        N.append(len(d["sim"][0])/2)
        err.append(d["RMS"])

    err = np.array(err)


    ax = fig.add_subplot(nr,nc,1); axes.append(ax)
    ax.plot(x, rho,'.',ms=1)
    ax.plot(x, r_rho,'-',lw=0.5)
    ax.set_ylabel(r"$\rho$")
    ax.set_xticklabels([])

    l2 = err[:,0]
    z = np.polyfit(np.log10(N),np.log10(l2),1)
    pz = np.poly1d(z)
    ax = fig.add_subplot(nr,nc,2)
    ax.plot(np.log10(N), np.log10(l2),'o',ms=5)
    ax.plot(np.log10(N), pz(np.log10(N)), "r--", alpha=0.5, label="O(%0.2g)"%z[0])
    ax.legend(fontsize=8, frameon=False)
    ax.set_xticklabels([])
    ax.set_ylabel(r"$\log \left(L_2\left(\rho\right)\right)$")

    ax = fig.add_subplot(nr,nc,3); axes.append(ax)
    ax.plot(x, u,'.',ms=1)
    ax.plot(x, r_u,'-',lw=0.5)
    ax.set_ylabel(r"$u$")
    ax.set_xticklabels([])

    l2 = err[:,1]
    z = np.polyfit(np.log10(N),np.log10(l2),1)
    pz = np.poly1d(z)
    ax = fig.add_subplot(nr,nc,4)
    ax.plot(np.log10(N), np.log10(l2),'o',ms=5)
    ax.plot(np.log10(N), pz(np.log10(N)), "r--", alpha=0.5, label="O(%0.2g)"%z[0])
    ax.legend(fontsize=8, frameon=False)
    ax.set_xticklabels([])
    ax.set_ylabel(r"$\log \left(L_2\left(u\right)\right)$")

    ax = fig.add_subplot(nr,nc,5); axes.append(ax)
    ax.plot(x, T,'.',ms=1)
    ax.plot(x, r_T,'-',lw=0.5)
    ax.set_ylabel(r"$T$")
    ax.set_xlabel(r"$x$")

    l2 = err[:,2]
    z = np.polyfit(np.log10(N),np.log10(l2),1)
    pz = np.poly1d(z)
    ax = fig.add_subplot(nr,nc,6)
    ax.plot(np.log10(N), np.log10(l2),'o',ms=5)
    ax.plot(np.log10(N), pz(np.log10(N)), "r--", alpha=0.5, label="O(%0.2g)"%z[0])
    ax.legend(fontsize=8, frameon=False)
    ax.set_ylabel(r"$\log \left(L_2\left(T\right)\right)$")
    ax.set_xlabel(r"$\log \left(N_\mathrm{cells}\right)$")

    for ax in axes:
        ax.set_xlim(-1,1)


    fig.tight_layout()

    fig.savefig("plot.pdf", dpi=300)
    plt.close(fig)

    return False
        

if __name__ == "__main__":
    sys.exit(check())