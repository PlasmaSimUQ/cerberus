
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
import matplotlib.ticker as ticker

import subprocess, shlex

def get_name(levels):
    return "Soliton_n%i.plt"%levels

def run(d):

    inputs = "Soliton.inputs"
    EXE = "../../local/MFP.1d.gnu.MPI.ex"

    ncells = d["n"]
    name = get_name(ncells)

    d["name"] = name

    if 1:
        run = "rm -r {base_name}*; mpirun -n 1 {exe} {input} amr.n_cell={n} amr.plot_file={base_name} 2>&1 | tee {base_name}.log".format(n=ncells, exe=EXE, base_name=name, input=inputs)
        subprocess.call(run, shell=True)

    return

def RMS(a, b):
    return np.sqrt(np.sum((a-b)**2))/float(len(a))

def get_data(f):

    rh5 = ReadBoxLib(f)
    t = rh5.time
        
    rh5 = ReadBoxLib(f, max_level=-1)

    name = "field"

    xc, By = rh5.get("y_B-%s"%name)
    xn, Dz = rh5.get("z_D-%s"%name, grid='node')

    rh5.close()

    return t, xc[0], By, Dz


def check():

    dom_len = 20

    data = [
        {"n":16*dom_len},
        {"n":32*dom_len},
        {"n":64*dom_len},
        {"n":128*dom_len},
    ]

    for d in data:

        run(d)


        # get a list of all the files in this directory
        files = get_files('.', include=[d["name"], '.plt'], exclude=["temp", "log"], get_all=True)

        t, x, By, Dz = get_data(files[-1])

        d["sim"] = [x, By, Dz]
            

        # =============================================================================
        # analytical solution
        # =============================================================================

        c = 1.0

        By_analytical = 0.5*(np.exp(-(-c*t + x)**2) + np.exp(-(c*t + x)**2))
        Dz_analytical = 0.5*c*(-np.exp(-(-c*t + x)**2) + np.exp(-(c*t + x)**2))

        d["ref"] = [x, By_analytical, Dz_analytical]

        # ERROR
        d["RMS"] = [RMS(By_analytical, By), RMS(Dz_analytical, Dz)]


    # =============================================================================
    # plot
    # =============================================================================

    plt.rc("font", family="serif")
    plt.rc("font", size=8)
    plt.rc("mathtext", fontset="cm")
    # matplotlib.rc('text', usetex = True)
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    error = False


    axes = []
    fig = plt.figure(figsize=(5,2.5))

    nr = 2
    nc = 2

    x, By, Dz = data[0]["sim"]
    x, r_By, r_Dz = data[0]["ref"]

    N = []
    err = []

    for d in data:
        N.append(len(d["sim"][0])/dom_len)
        err.append(d["RMS"])

    err = np.array(err)


    ax = fig.add_subplot(nr,nc,1); axes.append(ax)
    ax.plot(x, By,'.',ms=1)
    ax.plot(x, r_By,'-',lw=0.5)
    ax.set_ylabel(r"$B_y$")
    ax.set_xticklabels([])

    l2 = err[:,0]
    z = np.polyfit(np.log10(N),np.log10(l2),1)
    pz = np.poly1d(z)
    ax = fig.add_subplot(nr,nc,2)
    ax.plot(np.log10(N), np.log10(l2),'o',ms=5)
    ax.plot(np.log10(N), pz(np.log10(N)), "r--", alpha=0.5, label="O(%0.2g)"%z[0])
    ax.legend(fontsize=8, frameon=False)
    ax.set_xticklabels([])
    ax.set_ylabel(r"$\log \left(L_2\left(B_y\right)\right)$")

    if z[0] > 3:
        error = True

    ax = fig.add_subplot(nr,nc,3); axes.append(ax)
    ax.plot(x, Dz,'.',ms=1)
    ax.plot(x, r_Dz,'-',lw=0.5)
    ax.set_ylabel(r"$D_z$")
    ax.set_xlabel(r"$x$")

    l2 = err[:,1]
    z = np.polyfit(np.log10(N),np.log10(l2),1)
    pz = np.poly1d(z)
    ax = fig.add_subplot(nr,nc,4)
    ax.plot(np.log10(N), np.log10(l2),'o',ms=5)
    ax.plot(np.log10(N), pz(np.log10(N)), "r--", alpha=0.5, label="O(%0.2g)"%z[0])
    ax.legend(fontsize=8, frameon=False)
    ax.set_ylabel(r"$\log \left(L_2\left(D_z\right)\right)$")
    ax.set_xlabel(r"$\log \left(N_\mathrm{cells}\right)$")

    if z[0] > 3:
        error = True

    for ax in axes:
        ax.set_xlim(-10,10)


    fig.tight_layout()

    fig.savefig("plot.pdf", dpi=300)
    plt.close(fig)

    return error
        

if __name__ == "__main__":
    sys.exit(check())