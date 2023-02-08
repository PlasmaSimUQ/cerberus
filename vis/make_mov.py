# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:42:38 2016

@author: daryl
"""

import os, shutil
import matplotlib as mpl

mpl.use('agg')

from glob import glob
import h5py
from get_boxlib import ReadBoxLib, get_files, get_time
import pylab as plt
from PIL import Image, ImageOps
from multiprocessing import Pool
import subprocess, shlex
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FixedFormatter
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from shifted_colormap import shiftedColorMap
import types
import copy
import gc
import hashlib
import shelve

#==============================================================================
#
#==============================================================================

#plt.rc('text', usetex=True)
#plt.rc("font", family="serif")
#plt.rc("font", size=10)

verbosity = 3

#==============================================================================
#
#==============================================================================

def stringify(d, exclude=[]):
    s = []
    for key, value in d.items():
        if key in exclude:
            continue
        if type(value) == dict:
            s += stringify(value, exclude)
        elif type(value) == list:
            for i, v in enumerate(value):
                s += stringify({str(i):v}, exclude)
        elif callable(value):
            s.append(value.__name__)
        else:
            s.append(repr(value))

    return s

#==============================================================================
#
#==============================================================================

def distance(a,b):
    return np.sqrt((a.x - b.x)**2 + (a.y - b.y)**2)

def is_between(a,c,b, epsilon=1e-12):
    return -epsilon < (distance(a, c) + distance(c, b) - distance(a, b)) < epsilon

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def harvest_data_repeater(frame, limit=10):
    """
    keep trying
    """
    
    if frame["q"]["cores"] < 2:
        data = harvest_data(frame)
        return data

    for i in range(limit):
        try:
            data = harvest_data(frame)
            return data
        except:
            continue

    print("failed to harvest ",frame["h5"])

    return {"valid":False}


def harvest_data(frame):
    """
    return the limits of the data we are wanting to plot
    """

    h5Name = frame["h5"]
    q = frame["q"]
    get = q["get"]
    data = {}
    
    print("harvest : ",h5Name)

    if type(q["xy_limits"]) == types.FunctionType:
        t = get_time(h5Name)
        window = q["xy_limits"](t)
    else:
        window = q["xy_limits"]
    
    rc = ReadBoxLib(h5Name, max_level=q["level"], limits=window)
    
    t = rc.time

    for s in get:

        # determine the name of the file we want to save to
        output_path = s["output_path"]
        h5_save = os.path.split(h5Name)[1]+"-%s-level=%i.h5"%(s["tag"], q["level"])
        h5_save = os.path.join(output_path,h5_save)

        if frame["names_only"]:
            data[s] = {"h5":h5_save}
            continue

        # check if that file exists
        if (not os.path.isfile(h5_save)) or q["force_data"]:

            save_data = s["func"](rc, s)

            save_data["time"] = t
            save_data["window"] = window

            if "value" in save_data:
                limits = [save_data["value"].min(), save_data["value"].max()]
            else:
                limits = [0,0]

            h5 = h5py.File(h5_save,'w')
            for name, dat in save_data.items():
                h5.create_dataset(name, data=dat)
            h5.create_dataset("limits", data=limits)
            h5.close()

        else:
            h5 = h5py.File(h5_save,'r')
            limits = h5["limits"][()]
            h5.close()
            if verbosity > 2:
                print("harvest : ", h5_save, "already exists")

        data[s["tag"]] = {"h5":h5_save, "min":limits[0], "max":limits[1]}

    rc.close()
    
    if verbosity > 2:
        print("harvested: ", os.path.split(h5Name)[1])
        
    data["valid"] = True

    return data

def make_frame_repeater(frame, limit=10):
    """
    keep trying to make the frame
    """
    if frame["valid"]:
    
        if frame["q"]["cores"] < 2:
            make_frame(frame)
        else:
            success = False
            for i in range(limit):
                try:
                    make_frame(frame)
                    success = True
                    break
                except:
                    continue
            if not success:
                for i, s in enumerate(frame["q"]["get"]):
                    print("failed at attempt to make frame ",frame["id"]," for ",s["tag"])
    else:
        print("frame ",frame["id"]," invalid")

    frame.clear()
    gc.collect()

    return

def make_frame(frame):
    """
    make a single frame
    """

    q = frame["q"]

    output_path = q["frame_path"]
    fname = "%s-level=%i-%0.5d.png"%(q["name"], q["level"], frame["id"])
    fname = os.path.join(output_path,fname)
    frame["fname"] = fname

    # check if that file exists
    if not q["force_frames"]:
        if not (not os.path.isfile(fname)):
            if verbosity > 2:
                print("pre-existing: ",os.path.split(fname)[1])
            return
            

    # open and package up all of the h5 datasets that have been requested

    data = {}

    for s in q["get"]:
        tag = s["tag"]
        data[tag] = h5py.File(frame[tag]["h5"],'r')

    q["plot"](frame, data, fname)

    for s in q["get"]:
        data[s["tag"]].close()

    data.clear()

    # trim the image
    im=Image.open(fname)
    im_ = im.convert("RGB")
    im_ = ImageOps.invert(im_)
    bb = list(im_.getbbox())

    # make sure bb is divisible by 2
    for i in range(len(bb)):
        bb[i] = 2*int(bb[i]/2.0)

    # grow the bbox a little
    bb[0] -= 2
    bb[1] -= 2
    bb[2] += 2
    bb[3] += 2

    im2=im.crop(bb)
    im2.save(fname)

    if verbosity > 2:
        print("made :",os.path.split(fname)[1])

    im.close()
    im2.close()

    return True

def get_data(q, names_only=False):
    
    get_all = False
    if isinstance(q["time_span"], list):
        if not q["time_span"]:
            get_all = True

    names = get_files(q["files_dir"], include=q["file_include"], exclude=q["file_exclude"], times=q["time_span"], get_all=get_all)

    # make directories if they don't exist

    frames_path = os.path.join(q["mov_save"], q["name"])
    try:
        os.makedirs(frames_path)
    except OSError:
        if not os.path.isdir(frames_path):
            raise
    q["frame_path"] = frames_path

    for s in q["get"]:
        output_path = os.path.join(q["mov_save"], s["tag"])
        try:
            os.makedirs(output_path)
        except OSError:
            if not os.path.isdir(output_path):
                raise
        s["output_path"] = output_path
    
    


    q["data"] = []
    i = 0
    for name in names:
        if name == None:
            continue
        d = {}
        d["h5"] = name
        d["id"] = i
        d["names_only"] = names_only
        d["q"] = q

        q["data"].append(d)

        i += 1

    # go through all of the files and harvest the data
    DATA = []
    if q["cores"] == 1:
        for d in q["data"]:
            DATA.append(harvest_data_repeater(d))
    else:
        pool = Pool(processes=q["cores"])
        DATA = pool.map(harvest_data_repeater, q["data"])
        pool.close()
        pool.join()

    for i, data in enumerate(DATA):
        q["data"][i].update(data)

    return

def normalize(Q):

    limits = []

    # smooth the limits, but otherwise leave them alone
    for q in Q:
        if not isinstance(q["normalize"],dict):
            continue
        if "smooth" not in q["normalize"]:
            continue

        # get all of the limits first
        for s in q["get"]:
                
            vmin = []
            vmax = []
            ignore = []

            for i, frame in enumerate(q["data"]):
                if not frame["valid"]:
                    ignore.append(True)
                vmin.append(frame[s["tag"]]["min"])
                vmax.append(frame[s["tag"]]["max"])
                ignore.append(False)

            vmin = np.ma.masked_where(ignore, vmin)
            vmax = np.ma.masked_where(ignore, vmax)

            # apply a smoothing function
            vmin, vmax = q["normalize"]["smooth"](vmin, vmax)

            for i, frame in enumerate(q["data"]):
                if not frame["valid"]:
                    continue
                frame[s["tag"]]["min"] = vmin[i]
                frame[s["tag"]]["max"] = vmax[i]



    # get all the data limits and normalize within a folder
    for q in Q:
        lim = {}
        if q["normalize"] == "none":
            limits.append(None)
            continue
        if q["normalize"] in ["single", "all"]:
            for s in q["get"]:
                lim[s["tag"]] = [np.inf, -np.inf]

                for frame in q["data"]:
                    if not frame["valid"]:
                        continue
                    lim[s["tag"]][0] = min(lim[s["tag"]][0], frame[s["tag"]]["min"])
                    lim[s["tag"]][1] = max(lim[s["tag"]][1], frame[s["tag"]]["max"])

            limits.append(lim)

    # normalize across folders for each quantity
    for i, q in enumerate(Q):
        if q["normalize"] != "all":
            continue

        for s in q["get"]:
            for j, qq in enumerate(Q):
                if qq["normalize"] != "all":
                    continue

                try:
                    limits[i][s["tag"]][0] = min(limits[i][s["tag"]][0], limits[j][s["tag"]][0])
                    limits[i][s["tag"]][1] = max(limits[i][s["tag"]][1], limits[j][s["tag"]][1])
                except:
                    pass

    for i, q in enumerate(Q):
        if q["normalize"] not in ["all", "single"]:
            continue
        for s in q["get"]:
            for frame in q["data"]:
                if not frame["valid"]:
                    continue
                frame[s["tag"]]["min"] = limits[i][s["tag"]][0]
                frame[s["tag"]]["max"] = limits[i][s["tag"]][1]

    return

def get_streaklines(Q):

    name = "trajectory"

    # get all the streaklines
    for q in Q:

        for s in q["get"]:

            if "get_streak" not in s:
                continue
            
            X = []
            Y = []
            I = []

            for frame in q["data"]:
                if not frame["valid"]:
                    continue

                data = frame[s["tag"]]

                h5 = h5py.File(data["h5"],'r+')

                if name in h5:
                    if q["redo_streaks"]:
                        del h5[name]
                    else:
                        continue

                idata = h5["i"][()]
                rdata = h5["r"][()]
                x = rdata[:,0]
                y = rdata[:,1]
                idx = idata[:,0]
                cpu = idata[:,1]

                tag = cpu + idx*(cpu.max()+1)

                sorted_x = np.zeros((tag.max()+1))*np.nan
                sorted_y = np.zeros(sorted_x.shape)*np.nan

                sorted_x[tag] = x
                sorted_y[tag] = y
                
                X.append(sorted_x)
                Y.append(sorted_y)
                
                grp = h5.create_group(name)
                for i in range(len(X)):
                    grp.create_dataset("x_%i"%i, data=X[i])
                    grp.create_dataset("y_%i"%i, data=Y[i])
            
                grp.attrs["num"] = len(X)

                h5.close()

    return

def get_particle_trajectories(h5, xy_limits):
    px = []
    py = []
    pi = []

    grp = h5["trajectory"]
    n = grp.attrs["num"]

    max_len = 0
    for i in range(n):
        px.append(grp["x_%i"%i][()])
        py.append(grp["y_%i"%i][()])

        max_len = max(px[-1].size, max_len)
    
    px_ = np.nan*np.zeros((n, max_len))
    py_ = np.nan*np.zeros((n, max_len))

    for i in range(n):
        px_[i, 0:px[i].size] = px[i]
        py_[i, 0:py[i].size] = py[i]
    
    px = np.ma.masked_where(np.isnan(px_), px_)
    py = np.ma.masked_where(np.isnan(py_), py_)

    if n <= 1:
        return px, py, pi

    # mask out any portions of the streak that cross periodic boundaries
    pdx = np.diff(px, axis=0)
    pdy = np.diff(py, axis=0)

    step = np.zeros(px.shape)
    step[0:-1,:] = np.sqrt(pdx**2 + pdy**2)
    step[-1,:] = step[-2,:]
    dom_size = min((xy_limits[1][0] - xy_limits[0][0]), (xy_limits[1][1] - xy_limits[0][1]))
    mask = step > dom_size/2.0

    px = np.ma.masked_where(mask, px)
    py = np.ma.masked_where(mask, py)     

    return px, py, pi


def make_frames(q):


    if q["force_frames"]:
        delete_frames(q)


    # now make the frames
    if q["cores"] == 1:
        for frame in q["data"]:
            make_frame_repeater(frame)
    else:
        pool = Pool(processes=q["cores"], maxtasksperchild=2)
        pool.map(make_frame_repeater, q["data"])
        pool.close()
        pool.join()

    return


def make_movie(q):

    if q["only_frames"]:
        return

    path = q["frame_path"]
    all_files = listdir_fullpath(path)

    base_name = None
    for name in all_files:
        grab = True
        for tag in [".png", "level=%i"%q["level"]]:
            if tag not in name:
                grab = False
        if grab:
            base_name = name
            break

    if base_name == None:
        raise RuntimeError("no png files")

    base_name = base_name[0:-10]

    mov = base_name + ".mp4"

    prm = {"framerate":q["framerate"],
            "pattern":base_name + "-%05d.png",
            "output":mov
            }

    #call = "avconv -framerate %i -f image2 -i "%q["framerate"]
    call = "ffmpeg -framerate {framerate} -i {pattern} -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {output}".format(**prm)

    proc=subprocess.Popen(shlex.split(call))

    proc.communicate()

    proc=subprocess.Popen(shlex.split("mv %s %s"%(mov, q["mov_save"])))

    return

def delete_frames(q):


    output_path = q["frame_path"]
    fname = "%s-level=%i-*.png"%(q["name"], q["level"])
    fname = os.path.join(output_path,fname)
    for f in glob(fname):
        os.remove(f)

    return

def delete_h5(q):

    get_data(q, names_only=True)

    for data in q["data"]:
        for s in data:
            try:
                os.remove(data[s]["h5"])
            except:
                pass

    return


def make_all(Q):   
    
    # set default options
    for q in Q:
        items = {"normalize":"single",
                 }

        for item, value in items.items():
            if item not in q:
                q[item] = value

        exclude = [
        "framerate",
        "normalize",
        "cores",
        "force_frames",
        "only_frames",
        "redo_streaks",
        "dpi",
        "mov_save",
        "plot",
        "name",
    ]

    for q in Q:
        print(q["files_dir"])

        rep = stringify(q, exclude=exclude)
        rep = repr(sorted(rep)).encode('UTF8') 

        name = q["mov_save"]+"/"+hashlib.md5(rep).hexdigest() + ".shelf"

        q_needs_closing = False
        if os.path.isfile(name) and not q["force_data"]:
            q_ = shelve.open(name)

            for k, v in q_.items():
                if k not in q:
                    q[k] = v

            q_.close()
            print(" --> un-shelved")
        else:
            get_data(q)
            shelf = shelve.open(name)
            shelf.update(q)
            shelf.close()
            print(" --> shelved")

    

    print("normalizing...")
    normalize(Q)
    print("normalizing...done")

    print("getting particle trajectories...")
    get_streaklines(Q)
    print("getting particle trajectories...done")

    for q in Q:
        print(q["files_dir"])
        make_frames(q)
        make_movie(q)

    return


if __name__ == "__main__":

    """
    an example of how to use this script
    """

    def get_velocity_magnitude(ds):
        x, u = ds.get("x_vel-air")
        x, v = ds.get("y_vel-air", grid='node')

        return {"x":x, "value":np.sqrt(u**2 + v**2)}

    def get_alpha(ds):
        x, a = ds.get("alpha-air", grid='node')
        return {"x":x, "value":a}

    def get_vfrac(ds):
        x, v = ds.get("vfrac-air")
        return {"x":x, "value":v}

    def get_particles(ds):
        idat, rdat =  ds.get_particles('air')
        return {"i":idat, "r":rdat}

    def plot(frame, data, output_name):

        xc = data["vfrac-air"]["x"]
        v = data["vfrac-air"]["value"]
        
        xn = data["alpha-air"]["x"]
        a = data["alpha-air"]["value"]
        vmin = frame["alpha-air"]["min"]
        vmax = frame["alpha-air"]["max"]

        limits = frame["q"]["xy_limits"]

        # particles
        px, py, pi = get_particle_trajectories(data["particles-air"], limits)

        fig = plt.figure()
        ax = fig.add_subplot(111)

        yn, xn = np.meshgrid(xn[1], xn[0])
        yc, xc = np.meshgrid(xc[1], xc[0])

        ax.pcolormesh(xn, yn, a, vmin=vmin, vmax=vmax)

        ax.contour(xc, yc, v, levels=[0.5], colors=['k'])

        l = 10 #streak length
        ax.plot(px[-l::], py[-l::], "k-", lw=0.5, alpha=0.5)
        ax.plot(px[-1], py[-1], "ko", ms=0.5, alpha=0.5)

        ax.set_xlim(limits[0][0], limits[1][0])
        ax.set_ylim(limits[0][1], limits[1][1])

        ax.set_aspect(1)
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        fig.tight_layout()
        fig.savefig(output_name, dpi=300, bbox_inches="tight")
        plt.close(fig)

        return



    dt = 0.005

    Q = []

    q = {}
    q["files_dir"] = "../Exec/testing/EB-cavity"
    q["level"] = -1

    # all the data we need to retrieve
    q["get"] = [
        {"func":get_velocity_magnitude, "tag":"velocity-air"},
        {"func":get_alpha, "tag":"alpha-air"},
        {"func":get_vfrac, "tag":"vfrac-air"},
        {"func":get_particles, "tag":"particles-air", "get_streak":True},
    ]

    # how to make a frame
    q["plot"] = plot
    q["name"] = "test"
    
    ##
    q["framerate"] = 20
    q["mov_save"] = q["files_dir"] + "/mov"
    q["offset"] = [0.0, 0.0]
    q["xy_limits"] = [[0.0, 0.0], [1.0, 1.0]]
    q["file_include"] = ["plt"]
    q["file_exclude"] = ["chk"]
    q["cores"] = 8
    q["time_span"] = [] #np.arange(0,0.01,0.005).tolist()
    q["force_data"] = False
    q["force_frames"] = True
    q["only_frames"] = False
    q["redo_streaks"] = False
    q["dpi"] = 300

    q["normalize"] = "all"

    Q.append(q)

    make_all(Q)

    print("DONE")
