
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files, parse_header

import numpy as np
import pylab as plt
from matplotlib.image import NonUniformImage
from multiprocessing import Pool


#==============================================================================
# 
#==============================================================================

# plt_file = str(sys.argv[1])
plt_file = "TRMI"

component = "electron"
window = [[-0.4, 0.0], [0.75,1.0]]

# get a list of all the files in this directory
files = get_files('.', include=[plt_file], exclude=["chk","temp", ".png", "inputs"], times=[], tol=1e-4, get_all=True)

# get tracer particle data

def get_tracer_data(f):
    ds = ReadBoxLib(f)
    print(f," @ ",ds.time)
    
    idata, rdata = ds.get_particles(component)
    idx = idata[:,0]
    cpu = idata[:,1]
    x = rdata[:,0]
    y = rdata[:,1]
    
    sort_index = np.argsort(idx + cpu*idx.max())
    idx = idx[sort_index]
    cpu = cpu[sort_index]
    x = x[sort_index]
    y = y[sort_index]

    return x, y, idx

nproc = 1

if nproc > 1:
    pool = Pool(processes=nproc)
    DATA = pool.map(get_tracer_data, files)
    pool.close()
    pool.join()
else:
    DATA = []
    for f in files:
        DATA.append(get_tracer_data(f))

X = []
Y = []
I = []

for x,y,i in DATA:
    X.append(x)
    Y.append(y)
    I.append(i)

X = np.array(X)
Y = np.array(Y)
I = np.array(I)

# mask anything that goes periodic

pdx = np.diff(X, axis=0)
pdy = np.diff(Y, axis=0)

step = np.zeros(X.shape)
step[0:-1,:] = np.sqrt(pdx**2 + pdy**2)
step[-1,:] = step[-2,:]
dom_size = min((x.max() - x.min()), (y.max() - y.min()))
mask = step > dom_size/2.0

X = np.ma.masked_where(mask, X)
Y = np.ma.masked_where(mask, Y)
    
    
ds = ReadBoxLib(files[-1], max_level=-1, limits=window)

x, z = ds.get("rho-%s"%component)

x, y = x


# =============================================================================
# 
# =============================================================================


fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)


# plot the contour
im = NonUniformImage(ax, interpolation='bilinear', extent=(window[0][0], window[1][0], window[0][1], window[1][1]),
                     cmap="viridis")
im.set_data(x, y, z.T)
ax.images.append(im)
plt.colorbar(im, label=r"$\rho_\mathrm{%s}$"%component, orientation="horizontal")

# plot particle traces
plt.plot(X,Y,'w',lw=0.5, alpha=0.25)
plt.plot(X[-1,:],Y[-1,:],'ko', ms=0.5, alpha=0.5)

ax.set_xlim(window[0][0], window[1][0])
ax.set_ylim(window[0][1], window[1][1])
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect(1)

fig.tight_layout()

fig.savefig(plt_file+".png", dpi=300)
plt.close(fig)
    
    
print("DONE")
