
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_hdf5_data import ReadHDF5

import numpy as np
from scipy import fftpack
from scipy import signal
import pylab as plt
from matplotlib.image import NonUniformImage
from multiprocessing import Pool


#==============================================================================
# 
#==============================================================================

# plt_file = str(sys.argv[1])
plt_file = "TRMI"

window = [[-0.5, 1.5], [0,1]]

# get a list of all the files in this directory
files = ReadHDF5.get_files('.', include=[plt_file], exclude=["temp", ".png", "inputs"], times=[], tol=1e-4, get_all=True)

# get tracer particle data

rh5 = ReadHDF5(files[-1], max_level=-1, limits=window)

x, y, Dx = rh5.expression("{x_D-field}")
x, y, Dy = rh5.expression("{y_D-field}")
x, y, Dz = rh5.expression("{z_D-field}")

z = np.sqrt(Dx**2 + Dy**2 + Dz**2)

rh5.close()

# =============================================================================
# 
# =============================================================================

ni, nj = z.shape

n = 8

dx = x[1]-x[0]

ff = np.zeros(z.shape)*np.nan

o = int(n/2)

for i in range(o,ni-o):
    for j in range(o,nj-o):

        grab = z[i-o:i+o,j-o:j+o]

        f, wx = signal.welch(grab,dx,axis=0,nperseg=n)
        f, wy = signal.welch(grab,dx,axis=1,nperseg=n)

        wx = np.sum(wx,axis=1)
        wy = np.sum(wy,axis=0)

        w = wx+wy

        I = np.argmax(w)

        # print("w=",w)
        # print("I=",I)

        max_pwr_freq = f[I]

        # print("max of %g @ f = %g"%(w[I],f[I]))

        ff[i,j] = max_pwr_freq


# f = fftpack.fft2(z)

# f = np.log(np.abs(f))

# =============================================================================
# 
# =============================================================================


fig = plt.figure(figsize=(6,6))

###
ax = fig.add_subplot(211)

im = NonUniformImage(ax, interpolation='bilinear', extent=np.ravel(window),
                     cmap="viridis")
im.set_data(x, y, z.T)
ax.images.append(im)
plt.colorbar(im, label=r"$\left| \mathbf{D} \right|$", orientation="vertical")

ax.set_xlim(window[0][0], window[0][1])
ax.set_ylim(window[1][0], window[1][1])
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect(1)

###
ax = fig.add_subplot(212)

im = NonUniformImage(ax, interpolation='bilinear', extent=np.ravel(window),
                     cmap="viridis")
im.set_data(x, y, ff.T)
ax.images.append(im)
plt.colorbar(im, label=r"$\bar{f}$", orientation="vertical")

ax.set_xlim(window[0][0], window[0][1])
ax.set_ylim(window[1][0], window[1][1])
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect(1)

fig.tight_layout()

fig.savefig(plt_file+"_freq.png", dpi=300)
plt.close(fig)
    
    
print("DONE")
