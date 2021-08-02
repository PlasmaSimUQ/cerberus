# -*- coding: utf-8 -*-
"""
read AMReX data files
"""

import numpy as np
import os
import re
import random, string
import json
import copy
import pprint
from collections import namedtuple
import  tarfile
import codecs

from scipy.interpolate import RectBivariateSpline, UnivariateSpline


# ==============================================================================
# utility
# ==============================================================================

get_float = r"[-+]?\d*\.?\d+[eE]?[-+]?\d*"
get_int = r"[0-9]+"


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def path_join(pieces, sep="/"):
    path = ""
    for i, p in enumerate(pieces):
        path += p
        if i < len(pieces) -1:
            path += sep
    
    return path

def parse_MFP(name, data):
    # grab extra parameters

    fid = open(os.path.join(name, "info"))
    data.update(json.load(fid))
    fid.close()
    

def parse_FAB_header(name, dim):

    try:
        fid = open(name+"_H")
    except:
        tar_name, fab_header = os.path.split(name)
        tar = tarfile.open(tar_name+".tar")
        base_path, fab_header_path = os.path.split(tar_name)
        fab_header = path_join([fab_header_path,fab_header])
        fid = tar.extractfile(fab_header+"_H")
        fid = codecs.getreader("utf-8")(fid)

    data = {}

    data["version"] = int(fid.readline().rstrip())

    data["how"] = int(fid.readline().rstrip())

    data["ncomp"] = int(fid.readline().rstrip())

    numbers = re.findall(get_int, fid.readline().rstrip())
    nghost = np.array([int(x) for x in numbers])
    if len(nghost) == 1:
        data["nghost"] = np.ones((dim), dtype=int)*nghost[0]

    # box array data
    numbers = re.findall(get_int, fid.readline().rstrip())
    data["maxbox"] = int(numbers[0])

    data["boxes"] = []

    for i in range(data["maxbox"]):
        numbers = re.findall(get_int, fid.readline().rstrip())
        numbers = [int(x) for x in numbers]
        domains = np.reshape(numbers,(3, dim))


        data["boxes"].append({"lo":domains[0], "hi":domains[1], "type":domains[2]})

    fid.readline() # skip line

    # fab on disk
    n_fab_on_disk = int(fid.readline().rstrip())

    box_map = []

    for i in range(n_fab_on_disk):
        s = fid.readline().rstrip()
        words = s.split()
        path, back = os.path.split(name)
        data["boxes"][i]["fab_path"] = os.path.join(path, words[1])
        data["boxes"][i]["fab_offset"] = int(words[2])

    
    fid.readline() # skip line

    # data min/max
    numbers = re.findall(get_int, fid.readline().rstrip())
    for i in range(int(numbers[0])):
        min_val = re.findall(get_float, fid.readline().rstrip())
        data["boxes"][i]["min"] = np.array([float(x) for x in min_val])

    fid.readline() # skip line

    numbers = re.findall(get_int, fid.readline().rstrip())
    for i in range(int(numbers[0])):
        min_val = re.findall(get_float, fid.readline().rstrip())
        data["boxes"][i]["max"] = np.array([float(x) for x in min_val])
    
    fid.close()

    return data

def parse_plot_header(name, fid, data):

    # data names

    num_items = int(fid.readline().rstrip())

    names = []
    
    for i in range(num_items):
        names.append(fid.readline().rstrip())
    data["names"] = {"all":names}

    # dimensionality
    data["dim"] = int(fid.readline().rstrip())

    # time 
    data["time"] = float(fid.readline().rstrip())

    # number of levels
    data["n_levels"] = int(fid.readline().rstrip()) + 1

    data["levels"] = []
    for i in range(data["n_levels"]):
        data["levels"].append({})

    # domain lo
    dom_lo = fid.readline().rstrip().split(" ")

    #domain hi
    dom_hi = fid.readline().rstrip().split(" ")

    for i in range(data["dim"]):
        dom_lo[i] = float(dom_lo[i])
        dom_hi[i] = float(dom_hi[i])

    data["dom_lo"] = np.array(dom_lo)
    data["dom_hi"] = np.array(dom_hi)

    # refinement ratios
    refine = re.findall(get_int, fid.readline().rstrip())

    data["ref_ratio"] = np.array([])
    if refine:
        data["ref_ratio"] = np.array([int(x) for x in refine], dtype=int)


    # per level domain [level, (lo, hi, type), (x,y,z)]
    numbers = re.findall(get_int, fid.readline().rstrip())
    numbers = [int(x) for x in numbers]
    domains = np.reshape(numbers,(data["n_levels"], 3, data["dim"]))

    for i in range(data["n_levels"]):
        data["levels"][i]["domain_lo"] = domains[i,0]
        data["levels"][i]["domain_hi"] = domains[i,1]

    # level steps
    numbers = re.findall(get_int, fid.readline().rstrip())
    level_steps = [int(x) for x in numbers]

    for i in range(data["n_levels"]):
        data["levels"][i]["steps"] = level_steps[i]

    # level cell size
    for i in range(data["n_levels"]):
        numbers = re.findall(get_float, fid.readline().rstrip())
        data["levels"][i]["dx"] = np.array([float(x) for x in numbers])

    # coordinate system
    data["coord_sys"] = int(fid.readline().rstrip())

    # bwidth (??)
    data["bwidth"] = int(fid.readline().rstrip())

    for i in range(data["n_levels"]):
        level = data["levels"][i]

        numbers = re.findall(get_float, fid.readline().rstrip())

        level["level_id"] = int(numbers[0])
        level["num_grids"] = int(numbers[1])
        level["time"] = float(numbers[2])

        level["level_step"] = int(fid.readline().rstrip())

        # grid = [[lo],[hi]]
        level["grids"] = []
        for igrid in range(level["num_grids"]):
            grid = []
            for idim in range(data["dim"]):
                numbers = re.findall(get_float, fid.readline().rstrip())
                grid.append([float(x) for x in numbers])
            level["grids"].append(np.array(grid))

        level["data_loc"] = {"all":fid.readline().rstrip()}
        level["data"] = {"all":parse_FAB_header(os.path.join(name, level["data_loc"]["all"]), data["dim"])}

    return

def parse_checkpoint_header(name, fid, data): 

    # dimensionality
    data["dim"] = int(fid.readline().rstrip())
    dim = data["dim"]

    # time 
    data["time"] = float(fid.readline().rstrip())

    # number of levels
    data["max_level"] = int(fid.readline().rstrip())

    data["finest_level"] = int(fid.readline().rstrip())

    data["n_levels"] = data["max_level"] + 1

    data["levels"] = []
    for i in range(data["n_levels"]):
        data["levels"].append({})

    # skip some stuff as it is duplicated later
    for i in range(data["n_levels"]+1):
        fid.readline() # skip line

    # refinement ratios
    numbers = re.findall(get_int, fid.readline().rstrip())

    rr = np.reshape(numbers, (data["max_level"], data["dim"]))

    data["ref_ratio"] = rr[:,0].astype(int) # assume same refinement in all directions

    # time stepping info
    dt_level = re.findall(get_float, fid.readline().rstrip())
    dt_min = re.findall(get_float, fid.readline().rstrip())
    n_cycle = re.findall(get_float, fid.readline().rstrip())
    level_steps = re.findall(get_float, fid.readline().rstrip())
    level_count = re.findall(get_float, fid.readline().rstrip())
    
    
    for i in range(data["n_levels"]):
        level = data["levels"][i]
        level["dt"] = dt_level[i]
        level["dt_min"] = dt_min[i]
        level["n_cycle"] = n_cycle[i]
        level["level_steps"] = level_steps[i]
        level["level_count"] = level_count[i]

    data["names"] = {}

    for ilevel in range(data["n_levels"]):
        fid.readline() # skip level counter
        
        level = data["levels"][ilevel]

        # geometry
        coordsys = re.findall(get_float, fid.readline().rstrip())

        dx = np.zeros((dim))
        for i in range(dim):
            dx[i] = float(coordsys[1 + dim + i])

        level["dx"] = dx
        
            
        # [lo, hi, ilo, ihi, type, periodic]
        limits = re.findall(get_float, fid.readline().rstrip())
        limits = [float(x) for x in limits]
        probDomain = np.reshape(limits[0:dim*2], (2, dim), order='F')

        level["domain_lo"] = probDomain[0]
        level["domain_hi"] = probDomain[1]
        

        if ilevel == 0:
            data["dom_lo"] = level["domain_lo"]
            data["dom_hi"] = level["domain_hi"]

        domain = np.reshape(limits[dim*2:dim*5], (3, dim))
        level["lo"] = domain[0]
        level["hi"] = domain[1]
        level["periodic"] = domain[1]

        numbers = re.findall(get_int, fid.readline().rstrip())
        maxbox = int(numbers[0])

        for s in range(maxbox):
            fid.readline() # skip line

        # skip all of the per state data, this should be the same for all
        numbers = re.findall(get_int, fid.readline().rstrip())
        n_desc = int(numbers[0])

        level["data_loc"] = {}
        level["data"] = {}

        for n in range(n_desc):

            state_name = "state_%i"%n

            if n < data["num_state"]:
                state_name = data[state_name]["name"]

            fid.readline() # skip domain
            numbers = re.findall(get_int, fid.readline().rstrip())
            n_box = int(numbers[0])

            for s in range(n_box + 4):
                fid.readline() # skip line

            old = int(fid.readline().rstrip())

            loc = fid.readline().rstrip()
            level["data_loc"][state_name] = loc
            level["data"][state_name] = parse_FAB_header(os.path.join(name, loc), data["dim"])

            if old == 2:
                fid.readline() # just ignore old data

            # get the names of the data
            
            if (ilevel == 0) and (n < data["num_state"]):
                names = []
                dat = data["state_%i"%n]
                for i, cname in enumerate(dat["cons_names"]):
                    names.append(cname+"-"+dat["name"])
                data["names"][dat["name"]] = names
    



    return

def parse_particle_header(name, data, state_name):

    header_path = os.path.join(name, "Header")

    try:
        fid = open(header_path)
    except:
        tar = tarfile.open(name+".tar")
        base_path, particles_name = os.path.split(name)
#        header_path = os.path.join(particles_name, "Header")
        header_path = particles_name+"/Header" # use linux style path
        fid = tar.extractfile(header_path)
        fid = codecs.getreader("utf-8")(fid)

    version = fid.readline().rstrip()

    if "Version_Two_Dot_Zero_double" not in version:
        raise RuntimeError("particle header not supported")


    data["particles"] = {}
    dat = data["particles"]

    dat["dim"] = int(fid.readline().rstrip())
    dim = dat["dim"]
    dat["n_real_extra"] = int(fid.readline().rstrip())
    dat["n_real"] = dat["n_real_extra"] + dim

    dat["real_names"] = []
    for i in range(dim):
        dat["real_names"].append("position_%s"%["x","y","z"][i])
    for i in range(dat["n_real_extra"]):
        dat["real_names"].append(fid.readline().rstrip())

    dat["n_int_extra"] = int(fid.readline().rstrip())
    dat["n_int"] = dat["n_int_extra"] + 2

    dat["int_names"] = ["particle_id", "particle_cpu"]
    for i in range(dat["n_int_extra"]):
        dat["int_names"].append(fid.readline().rstrip())

    # skip checkpoint flag
    fid.readline() # skip line

    dat["n_particles"] = int(fid.readline().rstrip())

    #skip restart info
    fid.readline() # skip line
    
    num_levels = int(fid.readline().rstrip()) + 1

    grids_per_level = []
    for ilevel in range(num_levels):
        grids_per_level.append(int(fid.readline().rstrip()))

    n_grids = np.product(grids_per_level)

    for ilevel in range(num_levels):
        for igrid in range(grids_per_level[ilevel]):
            numbers = fid.readline().rstrip().split()
            ibox, npart, f_offset = [int(x) for x in numbers] #[file number, number of particles, offset]

            # if we actually have particles, then get where they are
            if npart > 0:

                pdata_name = os.path.join(name, "Level_%i"%ilevel, "DATA_%05i"%ibox)

                level_data = data["levels"][ilevel]["data"]

                if state_name in level_data:
                    level_data = level_data[state_name]
                else:
                    level_data = level_data["all"]
                
                box = level_data["boxes"][igrid]

                if "particles" not in box:
                    box["particles"] = []

                box["particles"].append({"state":state_name, "location":pdata_name, "number":npart, "offset":f_offset})

    return

def parse_header(name):

    data = {}

    parse_MFP(name, data)

    fid = open(os.path.join(name, "Header"))

    version = fid.readline().rstrip()

    if "HyperCLaw-V1.1" in version:
        parse_plot_header(name, fid, data)
    elif "CheckPointVersion_1.0" in version:
        parse_checkpoint_header(name, fid, data)
    else:
        raise RuntimeError("Header version not recognised")

    fid.close()
    # check for any particle data

    num_state = data["num_state"]
    for istate in range(num_state):
        state = data["state_%i"%istate]
        if "particle_data" in state:
            parse_particle_header(os.path.join(name, state["particle_data"]), data, state["name"])
        

    return data

def get_time(name):
    data = parse_header(name)
    return data["time"]

def get_files(folder, include=[], exclude=[], times=[], tol=1e-2, get_all=False):
    files = listdir_fullpath(folder)
    names = []
    for name in files:
        accept = True
        for exc in exclude:
            if exc in name:
                accept = False
                break
        for inc in include:
            if inc not in name:
                accept = False
                break

        if accept:
            names.append(name)

    names = sorted(names)

    if get_all or (times == []):
        return names

    file_times = []
    for idx, name in enumerate(names):
        try:
            data = parse_header(name)
        except:
            print("could not open " + name)
            pass

        file_times.append(data["time"])

    file_times = np.array(file_times)

    files = []

    # if times is an array
    try:
        times = times.tolist()
    except:
        pass

    if type(times) != list:
        times = [times]


    for time in times:
        if time < 0.0:
            files.append(names[int(time)])
            continue

        err = np.abs(file_times - time)

        fid = np.argmin(err)

        if err[fid] <= tol:
            files.append(names[fid])
        else:
            files.append(None)

    return files

def get_particle_files(folders):
    """
    given a list of folders retrieve the sub-folders that contain particle data
    """

    particles = {}


    for folder in folders:
        ppath = os.path.join(folder,"Particles")
        if not os.path.isdir(ppath):
            continue

        names = os.listdir(ppath)

        for name in names:
            path = os.path.join(ppath,name)
            if name in particles:
                particles[name].append(path)
            else:
                particles[name] = [path]


    return particles

# ==============================================================================
#
# ==============================================================================

class ReadBoxLib:

    def __init__(self, dataName, max_level=-1, limits=[]):

        self.dataName = dataName

        self.data = parse_header(dataName)

        # pprint.pprint(self.data["levels"][0])

        if limits:
            self.limits = copy.deepcopy(limits)
        else:
            self.limits = [self.data["dom_lo"], self.data["dom_hi"]]        

        n_levels = self.data["n_levels"]
        self.dim = self.data["dim"]

        if max_level < 0: 
            self.max_level = n_levels + max_level
        else:
            self.max_level = min(max_level, n_levels)


        if self.data:
            self.time = self.data["time"]
            self.names = []

            self.data["hydro_names"] = []
            self.data["hydro_mass"] = []
            self.data["hydro_charge"] = []
            self.data["hydro_gamma"] = []

            for i in range(self.data["num_state"]):
                istate = self.data["state_%i"%i]
                self.names.append(istate["name"])

                if istate["type"] == "hydro":
                    self.data["hydro_names"].append(istate["name"])
                    self.data["hydro_mass"].append(istate["mass"])
                    self.data["hydro_charge"].append(istate["charge"])
                    self.data["hydro_gamma"].append(istate["gamma"])

        self.flat_data = {}
        self.retrieval_sites = []
        self.refinement = []
        self.real_boxes = []

        return

    def retrieve(self):
        """
        generate a list of where to get data from
        """

        get_this = []

        dom_lo = self.data["dom_lo"]

        lo = self.limits[0][0:self.dim]
        hi = self.limits[1][0:self.dim]

        # go through all of the boxes and 
        # see if the box intersect the bounds

        for ilevel, level in enumerate(self.data["levels"]):

            # if ilevel > self.max_level:
            #     break

            dx = level["dx"]

            ilo = np.floor((lo - dom_lo + dx/2.0)/dx).astype(int)
            ihi = np.floor((hi - dom_lo - dx/2.0)/dx).astype(int)

            for istate, state in level["data"].items():

                # exclude anything that isn't part of a state definition
                if istate not in self.names and istate != "all":
                    continue

                for ibox, box in enumerate(level["data"][istate]["boxes"]):

                    box_lo = box["lo"]
                    box_hi = box["hi"]

                    box_real = [box_lo*dx + dom_lo, (box_hi+1)*dx + dom_lo]

                    intersect = np.all((box_lo <= ihi) & (box_lo >= ilo)) | np.all((box_hi <= ihi) & (box_hi >= ilo))

                    if intersect:
                        self.data["levels"][ilevel]["data"][istate]["boxes"][ibox]["fab_data"] = {} # initialise a place to put the data
                        ncomp = level["data"][istate]["ncomp"]
                        nghost = level["data"][istate]["nghost"]
                        get_this.append({"level_idx":ilevel, 
                            "box_idx":ibox, 
                            "ncomp":ncomp, 
                            "box_info":box, 
                            "nghost":nghost,
                            "state_idx":istate})
                        self.real_boxes.append(box_real)

        # now we have a list of where to get data

        self.retrieval_sites = get_this

        return

    def load_data(self, info, component):

        level_idx = info["level_idx"]

        if level_idx > self.max_level:
            return None

        box_idx = info["box_idx"]
        state_idx = info["state_idx"]

        fab_data = self.data["levels"][level_idx]["data"][state_idx]["boxes"][box_idx]["fab_data"]

        if component in fab_data:
            return fab_data[component]


        # get the component index
        if component in self.data["names"][state_idx]:
            comp_idx = self.data["names"][state_idx].index(component)
        else:
            return None
        

        binfo = info["box_info"]

        path = binfo["fab_path"]
        offset = binfo["fab_offset"]

        nghost = info["nghost"]

        lo = binfo["lo"] - nghost
        hi = binfo["hi"] + nghost
        shape = hi - lo + 1

        n_dat = np.product(shape)

        n_bytes = np.dtype("float64").itemsize

        try:
            f = open(path, "rb")
        except:
            path_pieces = path.split(os.path.sep)
            if path[0] == os.path.sep:
                path_pieces.insert(0,os.path.sep)
            tar_name = os.path.join(*(path_pieces[0:-1]))
            tar = tarfile.open(tar_name+".tar")
            data_path = path_join(path_pieces[-2::])
            f = tar.extractfile(data_path)

        f.seek(offset)
        f.readline()  # skip the first line
        f.seek(n_dat*comp_idx*n_bytes,1) # offset for component number
        buf = f.read(n_dat*n_bytes)
        arr = np.frombuffer(buf, "float64", n_dat)
        arr = arr.reshape(shape, order="F")

        slicer = []
        for n in nghost:
            if n == 0:
                slicer.append(slice(None))
            else:
                slicer.append(slice(n,-n))
        slicer = tuple(slicer)

        arr = arr[slicer]

        level_idx = info["level_idx"]
        box_idx = info["box_idx"]

        self.data["levels"][level_idx]["data"][state_idx]["boxes"][box_idx]["fab_data"][component] = arr

        return arr

    def refine_array(self, a, r):
        """
        refine an array by a factor of r
        uses linear slopes in sx and sy
        """

        o = a.repeat(r, axis=0)
        if self.dim == 1:
            ni = a.size
        elif self.dim == 2:
            o = o.repeat(r, axis=1)
            ni, nj = a.shape
            dj = 1/float(nj)
            jj = np.linspace(dj/2.0, 1-dj/2.0, nj)

        di = 1/float(ni)
        ii = np.linspace(di/2.0, 1-di/2.0, ni)

        if self.dim == 1:
            f = UnivariateSpline(ii, a, k=3, s=0)
        elif self.dim == 2:
            f = RectBivariateSpline(ii, jj, a, kx=3, ky=3, s=0)
            nj *= r
            dj = 1/float(nj)
            jj = np.linspace(dj/2.0, 1-dj/2.0, nj)

        ni *= r
        di = 1/float(ni)
        ii = np.linspace(di/2.0, 1-di/2.0, ni)

        if self.dim == 1:
            o = f(ii)
        elif self.dim == 2:
            o = f(ii, jj)

        return o
    
    def get_flat_limits(self):
        # cell spacing at the maximum level
        dx = self.data["levels"][self.max_level]["dx"]

        # indexing into finest level that matches the passed in limits
        dom_lo = self.data["dom_lo"]
        lo = self.limits[0][0:self.dim]
        hi = self.limits[1][0:self.dim]
        
        ix_lo = np.floor((lo - dom_lo + dx/2.0)/dx).astype(int)
        ix_hi = np.floor((hi - dom_lo + dx/2.0)/dx).astype(int)
        nx = ix_hi - ix_lo

        xc = []
        xn = []
        for i in range(self.dim):
            start = (0.5 + ix_lo[i])*dx[i] + dom_lo[i]
            stop  = ix_hi[i]*dx[i] + dom_lo[i]
            xc.append(np.arange(start, stop, dx[i]))

            start = ix_lo[i]*dx[i] + dom_lo[i]
            stop  = (0.5+ix_hi[i])*dx[i] + dom_lo[i]
            xn.append(np.arange(start, stop, dx[i]))

        return ix_lo, nx, xc, xn

    def get_flat(self, component, get_refinement=False):
        """
        return data on a uniform grid at cell spacing according to max_level
        """

        ix_lo, nx, xc, xn = self.get_flat_limits()

        if component in self.flat_data:
            return xc, self.flat_data[component]

        # make a repository for our final data
        fine = np.empty(nx)

        if get_refinement:
            refinement = np.zeros(nx)

        # collect the data

        ref_ratio = self.data["ref_ratio"].tolist()

        got_data = False

        dat = []
        for site in self.retrieval_sites:
            dat = self.load_data(site, component)

            if dat is None:
                continue
            else:
                got_data = True

            level = site["level_idx"]
            binfo = site["box_info"]

            # convert to the correct level indexing
            m = ref_ratio[level:self.max_level]
            m = int(np.prod(m))

            lo = binfo["lo"]*m - ix_lo
            hi = (binfo["hi"] + 1)*m - ix_lo


            if m > 1:
                dat = self.refine_array(dat, m)

            slicer_fine = []
            slicer_dat = []
            for i in range(self.dim):

                lo_i = max(lo[i], 0)
                hi_i = min(hi[i], nx[i])
                di = hi_i - lo_i
                i0 = lo_i
                i1 = i0 + di



                slicer_fine.append(slice(i0,i1))
                slicer_dat.append(slice(lo_i - lo[i], lo_i - lo[i] + di))

            fine[tuple(slicer_fine)] = dat[tuple(slicer_dat)]

            if get_refinement:
                refinement[tuple(slicer_fine)] = level

        if not got_data:
            raise RuntimeError("it appears that component '%s' doesn't exist, try one of %s"%(component, str(np.ravel(self.data["names"]))))

        self.flat_data[component] = fine

        if get_refinement:
            return xc, fine, refinement
        else:
            return xc, fine

    def mixture(self, alpha, name):
        names = self.data["hydro_names"]
        mass = self.data["hydro_mass"]
        charge = self.data["hydro_charge"]
        gamma = self.data["hydro_gamma"]
        
        sid = names.index(name)
        
        m0 = mass[sid][0]
        m1 = mass[sid][1]
        
        q0 = charge[sid][0]
        q1 = charge[sid][1]
        
        g0 = gamma[sid][0]
        g1 = gamma[sid][1]
        
        omass =   (m0*m1)/(m0*alpha + m1*(1.0-alpha))
        ocharge = (alpha*m0*q1 + (1.0-alpha)*m1*q0)/(m0*alpha + m1*(1.0-alpha))
        
        cp0 = g0/(m0*(g0-1.0))
        cp1 = g1/(m1*(g1-1.0))
        
        cv0 = 1.0/(m0*(g0-1.0))
        cv1 = 1.0/(m1*(g1-1.0))
        
        ogamma = ((1.0-alpha)*cp0 + alpha*cp1)/((1.0-alpha)*cv0 + alpha*cv1)
        
        self.flat_data["mass-"+name] = omass
        self.flat_data["charge-"+name] = ocharge
        self.flat_data["gamma-"+name] = ogamma
      
    def get(self, component, grid="cell", get_refinement=False):
        """
        get data corresponding to the names passed in
        """

        if not self.retrieval_sites:
            self.retrieve()
        
        #intercept special cases here
        for special in ["mass","charge","gamma"]:
            if special in component:

                # handle primitive vs conserved
                try:
                    x, alpha = self.get(component.replace(special, "alpha"))
                except:
                    x, trace = self.get(component.replace(special, "tracer"))
                    x, rho = self.get(component.replace(special, "density"))
                    alpha = trace/rho
                    
                name = component.replace(special,"")
                name = name.replace("-","")
                self.mixture(alpha, name)
                break

        flat = self.get_flat(component, get_refinement)
                
        if get_refinement:
            x, dat, ref = flat
        else:
            x, dat = flat

        if grid == "node":
            for i, x_ in enumerate(x):
                dx = x_[1] - x_[0]
                x[i] = np.linspace(x_[0] - dx/2.0, x_[-1]+dx/2.0, len(x_)+1)

        if get_refinement:
            return x, dat, ref
        else:
            return x, dat

    def load_particles(self, info, state_name):

        level_idx = info["level_idx"]
        box_idx = info["box_idx"]
        state_idx = info["state_idx"]

        if state_name != state_idx and state_idx != "all":
            return None

        box = self.data["levels"][level_idx]["data"][state_idx]["boxes"][box_idx]

        if "particles" not in box:
            return None

        part_data = box["particles"]

        n_real = self.data["particles"]["n_real"]
        n_int  = self.data["particles"]["n_int"]

        n_int_bytes = np.dtype("int32").itemsize
        n_real_bytes = np.dtype("float64").itemsize

        int_data = []
        real_data = []

        for i, grab in enumerate(part_data):

            if "data" in part_data:
                idata, rdata = grab["data"]
            else:
                npart = grab["number"]
                off = grab["offset"]
                path = grab["location"]
                state = grab["state"]

                if state_name != state:
                    continue

                try:
                    f = open(path, "rb")
                except:
                    path_pieces = path.split(os.path.sep)
                    if path[0] == os.path.sep:
                        path_pieces.insert(0,os.path.sep)
                    tar_name = os.path.join(*(path_pieces[0:-2]))
                    tar = tarfile.open(tar_name+".tar")
                    data_path = path_join(path_pieces[-3::])
                    f = tar.extractfile(data_path)

                # get the actual data
                f.seek(off)
                buf = f.read(n_int*npart*n_int*npart*n_int_bytes + n_real*npart*n_real_bytes)
                idata = np.frombuffer(buf, "int32", n_int*npart)
                idata = idata.reshape((npart, n_int), order="C")
                rdata = np.frombuffer(buf, "float64", n_real*npart, n_int*npart*n_int_bytes)
                rdata = rdata.reshape((npart, n_real), order="C")

                part_data[i]["data"] = (idata, rdata)

            int_data.append(idata)
            real_data.append(rdata)

        if len(int_data) > 0:
            int_data = np.concatenate(int_data)
            real_data = np.concatenate(real_data)

            return int_data, real_data
        else:
            return None

    def get_particles(self, state_name):

        if not self.retrieval_sites:
            self.retrieve()

        got_data = False

        int_data = []
        real_data = []

        for site in self.retrieval_sites:
            dat = self.load_particles(site, state_name)

            if dat is None:
                continue
            else:
                got_data = True

            int_data.append(dat[0])
            real_data.append(dat[1])

        if not got_data:
            raise RuntimeWarning("it appears that particles don't exist for '%s'"%(state_name))
        else:

            int_data = np.concatenate(int_data)
            real_data = np.concatenate(real_data)

            # remove entries that are not part of the requested interval
            
            ix_lo, nx, xc, xn = self.get_flat_limits()

            mask = (real_data[:,0] >= xn[0][0]) & (real_data[:,0] <= xn[0][-1])
            for i in range(1,self.dim):
                mask = mask & (real_data[:,i] >= xn[i][0]) & (real_data[:,i] <= xn[i][-1])

            int_data = int_data[mask,:]
            real_data = real_data[mask,:]

        return int_data, real_data

    def get_boxes(self):
        return self.real_boxes

    def close(self):
        self.data.clear()
        self.flat_data.clear()
        del self.retrieval_sites[:]


# ==============================================================================
#
# ==============================================================================

if __name__ == "__main__":
    import pylab as plt

    if 0:
        ds = ReadBoxLib("../Exec/testing/Viscous-Vortex/chk.Vortex00035")
        x, u = ds.get("nrg-air")

        plt.contourf(x[0], x[1], u)
        plt.gca().set_aspect(1)

    if 0:
        ds = ReadBoxLib("/home/uqdbond1/CODE/cerberus/Exec/testing/Riemann-MHD/Riemann.plt00095")
        x, u = ds.get("density-MHD")
        x, u = ds.get("y_B-MHD")

        plt.plot(x[0], u)

    if 1:
        f = "/home/uqdbond1/CODE/cerberus/Exec/testing/Isotope-RMI/IRMI.plt00307"

        ds = ReadBoxLib(f, max_level=-1, limits=[[-0.4, 0.0], [0.75,1.0]])

        idata, rdata = ds.get_particles("ion")
        
        x, u = ds.get("density-ion", grid="node")

        y, x = np.meshgrid(x[1], x[0])

        plt.pcolormesh(x, y, u)
        plt.plot(rdata[:,0], rdata[:,1], 'ko', ms=1)
        plt.gca().set_aspect(1)
    
        plt.show()
    
