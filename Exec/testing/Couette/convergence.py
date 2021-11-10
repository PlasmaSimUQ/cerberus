
import os, sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
from scipy.special import j1, y1
from numpy import exp, pi
from scipy.optimize import fsolve
import pylab as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import subprocess
import pickle


def make_inputs(template, name, options={}):
        tmp = open(template,'r')

        inputs_name =  name+".inputs"

        inp = open(inputs_name,"w")

        for line in tmp.readlines():
            try:
                line = line.format(**options)
            except:
                pass

            inp.write(line)

        
        tmp.close()
        inp.close()

        return inputs_name

def batch_run(doit=True):

    EXE="../../local/MFP.2d.gnu.MPI.EB.ex"

    resolutions = range(96,176+8,16)

    out = []

    for r in resolutions:

            name = "Couette_%i"%r

            inputs_name = make_inputs("Couette.inputs.template", name, {"nx":r, "ny":r, "name":name})

            if doit:
                if r < 32:
                    nproc=4
                else:
                    nproc=10
                run = "rm -r {base_name}.chk*; rm -r {base_name}.plt*; mpirun -n {nproc} {exe} {name} 2>&1 | tee {name}.log".format(exe=EXE, nproc=nproc, base_name=name, name=inputs_name)
                print(run)

                subprocess.call(run, shell=True)

            out.append({"r":r, "name":name})
    
    return out



def rootsearch(f,a,b,dx):
    """
    find if a root lies in the interval, but doesn't actually find it
    """
    x1 = a; f1 = f(a)
    x2 = a + dx; f2 = f(x2)
    while f1*f2 > 0.0:
        if x1 >= b:
            return None,None
        x1 = x2; f1 = f2
        x2 = x1 + dx; f2 = f(x2)
    return x1,x2

def roots(f, a, b, N, eps=1.0):
    """
    Returns a list of roots for the given function f
    The bounds used are shifted over the function by eps. 
    That is it begins with bounds x1 = x1 and x2 = x1 + eps
    """

    # first, generate a list of where we should find roots
    root_list = []
    while len(root_list) < N:
        x1,x2 = rootsearch(f,a,b,eps)

        x2 = min(x2, b)

        if x1 != None:
            a = x2
            root_list.append(0.5*(x1+x2))
        else:
            break

    # now solve for roots using scipy
    if root_list:

        root_list = fsolve(f, root_list)

    return root_list

def basis(n, R, R1):
    return j1(n*R)*y1(n*R1) - y1(n*R)*j1(n*R1)

def omega(R, R1, R2, omega_1, omega_2, nu, n, t):

        S0 = j1(n*R2)/(j1(n*R1)**2 - j1(n*R2)**2)*basis(n,R,R1)*(j1(n*R1)-j1(n*R2)*omega_1*R1/(omega_2*R2))

        S0 = pi*omega_2*R2*np.sum(S0)

        S1 =  j1(n*R2)/(j1(n*R1)**2 - j1(n*R2)**2)*basis(n,R,R1)*exp(-nu*n**2*t)*(omega_2*R2*j1(n*R1)-j1(n*R2)*omega_1*R1)

        S1 = -pi/R*np.sum(S1)

        return S0/R + S1

def check():

    plt.rc("font", family="serif")
    plt.rc("font", size=8)
    plt.rc("mathtext", fontset="cm")
    # matplotlib.rc('text', usetex = True)
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    # =============================================================================
    # roots for analytical solution
    # =============================================================================

    # D. J. Munoz, V. Springel, R. Marcus, M. Vogelsberger, L. Hernquist, 
    # Multidimensional, compressible viscous flow on a moving Voronoi mesh, 
    # Monthly Notices of the Royal Astronomical Society, 
    # Volume 428, Issue 1, 1 January 2013, Pages 254-279, 
    # https://doi.org/10.1093/mnras/sts015

    omega_1 = 0.5
    omega_2 = 0.1

    R1 = 1.0
    R2 = 2.5

    rho = 1.0
    p = 1.0
    mu = 0.005
    nu = mu/rho

    if os.path.isfile("roots.npy"):
        with open('roots.npy', 'rb') as f:
            n = np.load(f)
    else:
        n = roots(lambda x : basis(x, R2, R1), 1, 10000, 2000, 1)
        with open('roots.npy', 'wb') as f:
            np.save(f,n)

    #==============================================================================
    # Simulation results
    #==============================================================================

    plot_sim = False
    overwrite_pickle = True
    run_sims = True

    sim_results_name = "results.pickle"

    if overwrite_pickle or plot_sim or run_sims:

        # run a batch of simulations (see function for details)
        sims = batch_run(run_sims)

        for sim in sims:

            if plot_sim:
                fig = plt.figure(figsize=(5,3))

                gs = gridspec.GridSpec(nrows=2, ncols=2, width_ratios=[1,2.1], wspace=0.05, hspace=0.1)

                ax1 = fig.add_subplot(gs[0,0])
                ax2 = fig.add_subplot(gs[1,0])
                ax3 = fig.add_subplot(gs[:,1])

            # get a list of all the files in this directory
            files = get_files('.', include=['plt', sim["name"]], exclude=["temp"], get_all=True)

            # skip = int(len(files)/5)
            # files = [files[1]] + files[-1::-skip][::-1][1::]

            files = [files[-1]]

            L2_err = []

            for f in files:

                # print(f)

                data = ReadBoxLib(f)
                time = data.time
                    
                data = ReadBoxLib(f, max_level=-1)

                xn, u = data.get("x_vel-air", grid="node")
                xc, v = data.get("y_vel-air")
                xc, vf = data.get("vfrac-air")
                boxes = data.get_boxes()

                # u = np.ma.masked_where(vf < 0.99, u)
                # v = np.ma.masked_where(vf < 0.99, v)

                velocity = np.sqrt(u**2 + v**2)

                yc, xc = np.meshgrid(xc[1], xc[0])

                R = np.sqrt(xc**2 + yc**2)

                R = np.ravel(R)
                vel = np.ravel(velocity)

                # remove non-fluid entries
                mask = (R > R1) & (R < R2)
                R = R[mask]
                vel = vel[mask]

                # sort 
                I = np.argsort(R)
                R = R[I]
                vel = vel[I]

                # convert to angular velocity
                O = vel/R 

                # analytical solution
                # select = range(0, len(R), max(int(len(R)), int(len(R))/500))
                # aR = R[select]

                aR = R

                aO = np.zeros(aR.shape)
                for i, R_ in enumerate(aR):
                    aO[i] = omega(R_, R1, R2, omega_1, omega_2, nu, n, time)

                err = aO - O

                # L2 = np.sqrt(np.sum(err**2))/len(err)

                L2 = np.sqrt(np.sum(err**2)/float(err.size))

                L2_err.append(L2)

                # plotting
                if plot_sim:
                    ax1.plot(R, O,'.', ms=0.5, label=r"$t=%i$, L2=%0.1e"%(time, L2))
                    ax1.plot(aR, aO, 'k--', lw=0.5)
                    ax2.semilogy(aR, err, ".", ms=0.5, alpha=0.75)


            if plot_sim:
                # plot velocity field
                nx = len(xn[0])
                ny = len(xn[1])

                yn, xn = np.meshgrid(xn[1], xn[0])

                slicer = (slice(int(nx/2),nx), slice(int(ny/2),ny))

                cm = ax3.pcolormesh(xn[slicer], yn[slicer], velocity[slicer])

                # plot boxes
                grid = []
                for box in boxes:
                    sz = box[1] - box[0]
                    rect = Rectangle((box[0][0], box[0][1]), sz[0], sz[1])
                    grid.append(rect)

                pc = PatchCollection(grid, facecolor='none', alpha=1.0, edgecolor='w', linewidth=0.25)

                # plot quiver
                slicer2 = (slice(int(nx/2),nx,10), slice(int(ny/2),ny,10))
                ax3.quiver(xc[slicer2], yc[slicer2], u[slicer2], v[slicer2], pivot='mid')

                # plot boundaries
                
                cs = ax3.contour(xc[slicer], yc[slicer], vf[slicer], levels=[0.5], linewidths=[0.5], colors=['k'])


                for line in cs.allsegs[0]:

                    x = line[:,0].tolist()
                    y = line[:,1].tolist()

                    if (x[1] < x[-1]):
                        x += [2.6, 2.6, 0]
                        y += [0, 2.6, 2.6]
                    else:
                        x += [0]
                        y += [0]

                    plt.fill(x, y, 'w')
            

                # Add collection to axes
                ax3.add_collection(pc)

                ax3.set_aspect(1)
                ax3.set_xticks([])
                ax3.set_yticks([])
                ax3_divider = make_axes_locatable(ax3)
                # add an axes to the right of the main axes.
                cax = ax3_divider.append_axes("right", size="7%", pad="2%")
                cb = plt.colorbar(cm, cax=cax, label=r"$\left| \vec{u} \right|$")

                ax3.text(0.95, 0.95, r"$t=%0.2g$"%time, horizontalalignment='right', verticalalignment='top', transform=ax3.transAxes)

                ax1.set_xlim(1,2.5)
                ax1.set_ylabel(r"$\Omega$")
                ax1.set_xticklabels([])
                # ax1.legend()

                ax2.set_xlim(1,2.5)
                ax2.set_xlabel(r"$r$")
                ax2.set_ylabel(r"$\left|\hat{\Omega} - \Omega\right|$")

                # fig.tight_layout()

                fig.savefig(sim["name"]+".png",dpi=300, bbox_inches='tight')

            sim["L2"] = L2_err

        if overwrite_pickle:
            with open(sim_results_name, 'wb') as f:
                pickle.dump(sims, f, protocol=pickle.HIGHEST_PROTOCOL)
        
    with open(sim_results_name, 'rb') as f:
        sims = pickle.load(f)


    # plot the L2 error
    fig = plt.figure(figsize=(3,2))
    ax = fig.add_subplot(111)

    rr = []

    L2_max = []

    for sim in sims:
        r = sim["r"]
        L2 = sim["L2"]

        ax.plot(r*np.ones((len(L2))), L2, 'ok', ms=2.0)

        # ax.errorbar([r], [np.median(L2)], color='k', lw=0.75,xerr=0, yerr=np.max(L2) - np.median(L2), capthick=0.75, capsize=2)

        rr.append(r)
        # L2_max.append(np.max(L2))
        L2_max.append(L2[-1])

    z = np.polyfit(np.log10(rr),np.log10(L2_max),1)
    p = np.poly1d(z)


    ax.plot(rr, 10**p(np.log10(rr)), "r--", lw=0.75, alpha=0.75, label=r"$\frac{\mathrm{d} L_2}{\mathrm{d}N_c} = %0.3g$"%z[0])

    ax.legend()


    ax.set_xscale('log')
    ax.set_yscale('log')

    import matplotlib.ticker as mticker

    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())

    ax.tick_params(axis = "x", which = "minor", bottom = False, top = False)

    ax.set_xticks(rr)

    ax.set_xlabel(r"$N_c$")
    ax.set_ylabel(r"$L_2 = \sqrt{\frac{\sum_{i=1}^{N_c} \left|\hat{\Omega}_i - \Omega_i\right|^2}{N_c}}$")

    fig.tight_layout()

    fig.savefig("Couette_L2.pdf",dpi=300, bbox_inches='tight')



    # if np.max(L2_err) > 1e-3:
    #     return 1
    # else:
    #     return 0





if __name__ == "__main__":
    sys.exit(check())
