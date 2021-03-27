import pylab as plt
from mpl_toolkits import mplot3d
import matplotlib as mpl
import numpy as np
import h5py

def get_range(vmin, vmax):
    big = max(abs(vmin), abs(vmax))
    vmin = -big
    vmax = big
    return vmin, vmax

def visualise_face_array(fx, fy, fz):
    
    fig = plt.figure()
    
    axes = []
    
    ax = fig.add_subplot(1,3,1, projection='3d'); axes.append(ax)
    visualise_array_x(fx, ax)
    
    ax = fig.add_subplot(1,3,2, projection='3d'); axes.append(ax)
    visualise_array_y(fy, ax)
    
    ax = fig.add_subplot(1,3,3, projection='3d'); axes.append(ax)
    visualise_array_z(fz, ax)
    
    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_aspect(1)
        
        
    fig.tight_layout()
        
    return
    
    

def visualise_array_x(arr, ax):
    vis_arr = np.ma.masked_invalid(arr)
#    print "x:"
#    print vis_arr
    sx, sy, sz = vis_arr.shape
    x, y, z = np.mgrid[0:sx, -0.5:sy-0.5:(sy+1)*1j, -0.5:sz-0.5:(sz+1)*1j]
    vmin, vmax = get_range(vis_arr.min(), vis_arr.max())
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    for i in range(sx):
        my_col = mpl.cm.bwr(vis_arr[i,:,:])
        ax.plot_surface(x[i,:,:], y[i,:,:], z[i,:,:], norm=norm, facecolors = my_col, alpha=0.5)
        
    ax.set_title("x")
        
    return

def visualise_array_y(arr, ax):
    vis_arr = np.ma.masked_invalid(arr)
#    print "y:"
#    print vis_arr
    sx, sy, sz = vis_arr.shape
    x, y, z = np.mgrid[-0.5:sx-0.5:(sx+1)*1j, 0:sy, -0.5:sz-0.5:(sz+1)*1j]
    vmin, vmax = get_range(vis_arr.min(), vis_arr.max())
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    for i in range(sy):
        my_col = mpl.cm.bwr(vis_arr[:,i,:])
        ax.plot_surface(x[:,i,:], y[:,i,:], z[:,i,:], norm=norm, facecolors = my_col, alpha=0.5)
    ax.set_title("y")
    return

def visualise_array_z(arr, ax):
    vis_arr = np.ma.masked_invalid(arr)
#    print "z:"
#    print vis_arr
    sx, sy, sz = vis_arr.shape
    x, y, z = np.mgrid[-0.5:sx-0.5:(sx+1)*1j, -0.5:sy-0.5:(sy+1)*1j, 0:sz]
    vmin, vmax = get_range(vis_arr.min(), vis_arr.max())
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    for i in range(sz):
        my_col = mpl.cm.bwr(vis_arr[:,:,i])
        ax.plot_surface(x[:,:,i], y[:,:,i], z[:,:,i], norm=norm, facecolors = my_col, alpha=0.5)
    ax.set_title("z")
    return

# =============================================================================
# 
# =============================================================================

def multi_slice_viewer(volume, ax):

    ax.volume = volume
    data = np.ma.masked_invalid(volume)
    ax.index = volume.shape[-1] / 2
    dat = volume[:,:,ax.index].T
    im = ax.imshow(dat, origin="lower", aspect='equal',
                   vmin=data.min(), vmax=data.max())
    im.cmap.set_bad('magenta')
    
    hi = dat.shape
    
    # Minor ticks
    ax.set_xticks(np.arange(-.5, hi[1]+0.5, 1), minor=True);
    ax.set_yticks(np.arange(-.5, hi[0]+0.5, 1), minor=True);
    
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='w', linestyle='-', linewidth=0.5)
    
    ax.figure.canvas.mpl_connect('scroll_event', process_key)
    plt.colorbar(im)
    
    return ax

def process_key(event):
    fig = event.canvas.figure
    ax = fig.axes[0]
    if event.button == 'up':
        previous_slice(ax)
    elif event.button == 'down':
        next_slice(ax)
    fig.canvas.draw()

def previous_slice(ax):
    """Go to the previous slice."""
    volume = ax.volume
    ax.index = (ax.index - 1) % volume.shape[-1]  # wrap around using %
    ax.images[0].set_array(volume[:,:,ax.index].T)
    print ax.index

def next_slice(ax):
    """Go to the next slice."""
    volume = ax.volume
    ax.index = (ax.index + 1) % volume.shape[-1]
    ax.images[0].set_array(volume[:,:,ax.index].T)
    print ax.index
    
def visualise_array_2d(fx):
    
    fig = plt.figure()
    
    ax = fig.add_subplot(1,1,1)
    multi_slice_viewer(fx, ax)
    
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect(1)
        
    fig.tight_layout()
        
    return

# =============================================================================
# 
# =============================================================================

def visualise_array_subset(arr, alo, ahi, lo, hi):
    
    arr = np.ma.masked_invalid(arr)
    
    sx, sy, sz = arr.shape
    x, y, z = np.mgrid[alo[0]:ahi[0]+1, alo[1]:ahi[1]+1, alo[2]:ahi[2]+1]
    
    
    fig = plt.figure()
    
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.scatter(x, y, z, s=100.0, c=arr.ravel())
    
    ####
    
    bounds = np.array([lo,hi], dtype=np.float)
    
    bounds[0,:] -= 0.5
    bounds[1,:] += 0.5
    
    pts = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pts.append((bounds[i,0], bounds[j,1], bounds[k,2]))
    
    box = [(0,4), (4,6), (6,2), (2,0), (1,5), (5,7), (7,3), (3,1),
           (0,1), (4,5), (6,7), (2,3)]
    
    for b in box:
        line = np.array([pts[b[0]],pts[b[1]]])
        ax.plot(line[:,0], line[:,1], line[:,2],'r')
    
    #ax.plot()
    
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect(1)
        
    fig.tight_layout()
    
    
    return

# =============================================================================
# 
# =============================================================================
    

def visualise_array_3d(arr, alo, ahi):
    
    arr = np.ma.masked_invalid(arr)
    
    sx, sy, sz = arr.shape
    x, y, z = np.mgrid[alo[0]:ahi[0]+1, alo[1]:ahi[1]+1, alo[2]:ahi[2]+1]
    
    
    fig = plt.figure()
    
    ax = fig.add_subplot(1,1,1, projection='3d')
    p = ax.scatter(x, y, z, s=100.0, c=arr.ravel())
    fig.colorbar(p)

# =============================================================================
# 
# =============================================================================
    
def save_state(cd, cx, cy, cz, clo, chi, cslo, cshi, fx, fy, fz, fc, flo, fhi, fslo, fshi, dx, name):
    
    h5 = h5py.File(name+".hdf5","w")
    
    h5.create_dataset("cd", data=cd)
    h5.create_dataset("cx", data=cx)
    h5.create_dataset("cy", data=cy)
    h5.create_dataset("cz", data=cz)
    
    h5.create_dataset("clo", data=clo)
    h5.create_dataset("chi", data=chi)
    h5.create_dataset("cslo", data=cslo)
    h5.create_dataset("cshi", data=cshi)
    
    h5.create_dataset("fx", data=fx)
    h5.create_dataset("fy", data=fy)
    h5.create_dataset("fz", data=fz)
    h5.create_dataset("fc", data=fc)
    
    h5.create_dataset("flo", data=flo)
    h5.create_dataset("fhi", data=fhi)
    h5.create_dataset("fslo", data=flo)
    h5.create_dataset("fshi", data=fhi)
    
    h5.create_dataset("dx", data=dx)
    
    
    h5.close()
    
    return

def analyse_state(name):
    
    h5 = h5py.File(name+".hdf5","r")
    
    cd = np.ma.masked_invalid(h5["cd"][()])
    cx = np.ma.masked_invalid(h5["cx"][()])
    cy = np.ma.masked_invalid(h5["cy"][()])
    cz = np.ma.masked_invalid(h5["cz"][()])
    
    clo = h5["clo"][()]
    chi = h5["chi"][()]
    
    fx = np.ma.masked_invalid(h5["fx"][()])
    fy = np.ma.masked_invalid(h5["fy"][()])
    fz = np.ma.masked_invalid(h5["fz"][()])
    fc = np.ma.masked_invalid(h5["fc"][()])
    
    flo = h5["flo"][()]
    fhi = h5["fhi"][()]
    
    slo = h5["slo"][()]
    shi = h5["shi"][()]
    
    dx = h5["dx"][()]
    
    h5.close()
    
#    visualise_array_subset(fx, flo, fhi+np.array([1,0,0]), slo, shi+np.array([1,0,0]))
    
#    visualise_array_subset(fc, flo, fhi, slo, shi)
#    visualise_array_subset(cd, clo, chi, clo, chi)
    
    idx = 0
    
    div = (fx[1::,:,:,idx] - fx[0:-1,:,:,idx])/dx[0] + (fy[:,1::,:,idx] - fy[:,0:-1,:,idx])/dx[1] + (fz[:,:,1::,idx] - fz[:,:,0:-1,idx])/dx[2]
    
    #visualise_array_3d(fx[:,:,:,idx], flo, fhi+np.array([1,0,0]))
    visualise_array_3d(div, flo, fhi)
    #visualise_array_3d(fc, flo, fhi)
          
    err = np.log10(np.abs(fc - div))
    
    visualise_array_3d(err, flo, fhi)
    
    
    return
    
    
# =============================================================================
#     
# =============================================================================


if __name__ == "__main__":
    
    #analyse_state("initial")
    analyse_state("final")
    plt.show()