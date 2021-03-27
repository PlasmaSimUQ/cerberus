
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
    files = get_files('.', include=['.plt'], exclude=["temp"], get_all=True)

    f = files[-1]

    ds = ReadBoxLib(f)

    t = ds.time

    x, rho, ratio = ds.get('rho-MHD', get_refinement=True)
    x, By = ds.get('y_B-MHD')
    x, Bz = ds.get('z_B-MHD')

    x = x[0]

    sim_data = [rho, By, Bz]

    # =============================================================================
    # analytical solution
    # =============================================================================

    

    r_x = t*np.array([-2.500000, -1.474922, -0.990247, -0.631585, -0.631585, 
            -0.521395, -0.445268, 0.402052, 0.402052, 1.279598, 
            1.279598,  1.568067, 1.568067, 2.072332, 2.072332, 
            2.500000])
    r_rho =np.array([ 3.000000, 3.000000, 2.340949, 2.340949, 2.340949, 2.340949, 
            2.200167, 2.200167, 1.408739, 1.408739, 1.054703, 1.054703, 1.054703, 
            1.054703, 1.000000, 1.000000])
    r_By =[1.000000, 1.000000, 0.642777, 0.642777, 0.344252, 0.344252, 
            0.413199, 0.413199, 0.413199, 0.413199, 0.601050, 0.601050, 
            0.079386, 0.079386, 0.070737, 0.070737]
    r_Bz =[0.000000, 0.000000, 0.000000, 0.000000, 0.542820, 0.542820, 
            0.651535, 0.651535, 0.651535, 0.651535, 0.947741, 0.947741, 
            1.119452, 1.119452, 0.997495, 0.997495]

    ref_data = [r_rho, r_By, r_Bz]

    # =============================================================================
    # check
    # =============================================================================

    sample_x = [
        [-0.6, -0.58, -0.4, -0.2, -0.12, 0.12, 0.2, 0.5, 0.52, 0.82, 0.84],
        [-0.6, -0.58, -0.4, -0.27, -0.15, 0.5, 0.52, 0.6, 0.65, 0.82, 0.84],
        [-0.27, -0.15, 0.5, 0.52, 0.6, 0.65, 0.82, 0.84],
    ]

    sample_y = []
    sample_err  = []

    success = 0

    for i in range(len(sim_data)):
        sim = sim_data[i]
        sy = np.interp(sample_x[i], x, sim)
        sample_y.append(sy)

        ref = ref_data[i]
        ry = np.interp(sample_x[i], r_x, ref)

        err = ry - sy

        Linf = np.max(np.abs(err))
        sample_err.append(err)

        if (Linf > 0.05):
            success = 1

    # =============================================================================
    # plot
    # =============================================================================

    n_skip = 1
    x = x[::n_skip]
    rho = rho[::n_skip]
    By = By[::n_skip]
    Bz = Bz[::n_skip]
    ratio = ratio[::n_skip]

    plt.rc("font", family="serif")
    plt.rc("font", size=8)
    plt.rc("mathtext", fontset="cm")
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    axes = []
    fig = plt.figure(figsize=(5,4))

    nr = 3
    nc = 1

    ax = fig.add_subplot(nr,nc,1); axes.append(ax)
    ax.plot(x, rho,'-',lw=0.5)
    ax.plot(r_x, r_rho,'--',lw=0.5)
    ax.set_ylabel(r"$\rho$")
    ax.set_xticklabels([])
    
    ax = ax.twinx()
    ax.plot(x, ratio, 'k-.', lw=0.25)
    ax.set_yticks([0,2,4])
    ax.set_ylabel(r"level")
    

    ax = fig.add_subplot(nr,nc,2); axes.append(ax)
    ax.plot(x, By,'-',lw=0.5)
    ax.plot(r_x, r_By,'--',lw=0.5)
    ax.set_ylabel(r"$B_y$")
    ax.set_xticklabels([])

    ax = fig.add_subplot(nr,nc,3); axes.append(ax)
    ax.plot(x, Bz,'-',lw=0.5)
    ax.plot(r_x, r_Bz,'--',lw=0.5)
    ax.set_ylabel(r"$B_z$")
    ax.set_xlabel(r"$x$")


    
    for i, ax in enumerate(axes):
        ax.set_xlim(x[0], x[-1])

        # plot where the sample points are
        # ax.plot(sample_x[i], sample_y[i], 'ro', mfc='none')

    fig.tight_layout()

    fig.savefig("plot.pdf", dpi=300)
    plt.close(fig)

    return success
        

if __name__ == "__main__":
    sys.exit(check())