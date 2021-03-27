
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
from matplotlib.image import NonUniformImage
import matplotlib.ticker as ticker


#==============================================================================
# 
#==============================================================================

stuff = {}

# get a list of all the files in this directory
files = get_files('.', include=['plt'])

# get tracer particle data

def get_tracer_data(f):
    ds = ReadBoxLib(f)
    print(f," @ ",ds.time)
    
    idata, rdata = ds.get_particles('air')

    return idata, rdata

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

# due to the invalidation of particles we have a varying number of particles
# throughout the simulation. Therefore, we need to match up our particles over time
# in order to plot the streak lines

# find the unique ids
idx = set()

for pairs in DATA[0][0]:
    idx.add(tuple(pairs))

traces = {}
for i in idx:
    x = []
    y = []
    

    for d in DATA:
        ii = np.where((d[0][:,0] == i[0]) & (d[0][:,1] == i[1]))
        if len(ii[0]) > 0:
            ii = ii[0][0]
        else:
            continue
        x.append(d[1][ii,0])
        y.append(d[1][ii,1])

    traces[i] = {"x":x, "y":y}

    

# get data

f = files[-1]

data = ReadBoxLib(f)
t = data.time
    
data = ReadBoxLib(f, max_level=-1)

# x, rho = data.get("density-air", grid="node")
# x, mx = data.get("x_mom-air", grid="node")
# x, my = data.get("y_mom-air", grid="node")
xc, vf = data.get("vfrac-air", grid="cell")

# rho = np.ma.masked_where(vf<=1e-3, rho)

# u = mx/rho
# v = my/rho
# vel = np.sqrt(u**2 + v**2)

# y, x = np.meshgrid(x[1], x[0])
yc, xc = np.meshgrid(xc[1], xc[0])

for idx, trace in traces.items():

    x = trace["x"]

    if len(x) < 2:
        continue

    
    y = trace["y"]

    pdx = np.diff(x)
    pdy = np.diff(y)

    step = np.zeros(len(x))
    step[0:-1] = np.sqrt(pdx**2 + pdy**2)
    step[-1] = step[-2]
    dom_size = min((xc.max() - xc.min()), (yc.max() - yc.min()))
    mask = step > dom_size/2.0

    x = np.ma.masked_where(mask, x)
    y = np.ma.masked_where(mask, y)

    traces[idx] = {"x":x, "y":y}

# plot stuff

fig = plt.figure()
ax = fig.add_subplot(111)
# ax.pcolormesh(x, y, np.ma.masked_where(vf<=1e-3, vel))
ax.contour(xc, yc, vf, levels=[0.5], colors=['k'])

for k, t in traces.items():
    ax.plot(t["x"], t["y"],'k-',lw=0.5, alpha=1.0)
    if len(t["x"]) > 1:
        ax.plot([t["x"][-1]],[t["y"][-1]],'ko', ms=0.5, alpha=1.0)

ax.set_aspect(1)
fig.savefig("plot.png", dpi=300)
# plt.show()
    
    
print("DONE")
