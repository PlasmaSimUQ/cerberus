
import sys  # nopep8
cmd_folder = "../../../vis"  # nopep8
if cmd_folder not in sys.path:  # nopep8
    sys.path.insert(0, cmd_folder)

import numpy as np
import subprocess, shlex
import pylab as plt
import matplotlib.gridspec as gridspec

from get_boxlib import ReadBoxLib, get_files
from MMS_sources import hydro_var, hydro_src, hydro_const
from MMS_sources import field_var, field_src, field_const
from MMS_sources import mhd_var, mhd_src, mhd_const

def check():

    components = {}

    for i, var in enumerate(hydro_var):
        H = hydro_var[i]
        S = hydro_src[i]
        components["rho-fluid_%i"%i] = {"f":H["f_rho_%i"%i], "s":S["s_rho_%i"%i]}
        components["x_mom-fluid_%i"%i]   = {"f":H["f_mx_%i"%i],  "s":S["s_mx_%i"%i]}
        components["y_mom-fluid_%i"%i]   = {"f":H["f_my_%i"%i],  "s":S["s_my_%i"%i]}
        components["z_mom-fluid_%i"%i]   = {"f":H["f_mz_%i"%i],  "s":S["s_mz_%i"%i]}
        components["nrg-fluid_%i"%i]     = {"f":H["f_nrg_%i"%i], "s":S["s_nrg_%i"%i]}

    for d in ["x", "y", "z"]:
        for f in ["D", "B"]:
            components["%s_%s-field"%(d, f)] = {"f":field_var["f_%s%s"%(f,d)], "s":field_src["s_%s%s"%(f,d)]}

    components["phi-field"] = {"f":field_var["f_phi"], "s":field_src["s_phi"]}
    components["psi-field"] = {"f":field_var["f_psi"], "s":field_src["s_psi"]}

    components["rho-mhd"] = {"f":mhd_var["f_mhd_rho"], "s":mhd_src["s_mhd_rho"]}
    components["x_mom-mhd"]   = {"f":mhd_var["f_mhd_mx"],  "s":mhd_src["s_mhd_mx"]}
    components["y_mom-mhd"]   = {"f":mhd_var["f_mhd_my"],  "s":mhd_src["s_mhd_my"]}
    components["z_mom-mhd"]   = {"f":mhd_var["f_mhd_mz"],  "s":mhd_src["s_mhd_mz"]}
    components["nrg-mhd"]     = {"f":mhd_var["f_mhd_nrg"], "s":mhd_src["s_mhd_nrg"]}
    components["x_B-mhd"]   = {"f":mhd_var["f_mhd_Bx"], "s":mhd_src["s_mhd_Bx"]}
    components["y_B-mhd"]   = {"f":mhd_var["f_mhd_By"], "s":mhd_src["s_mhd_By"]}
    components["z_B-mhd"]   = {"f":mhd_var["f_mhd_Bz"], "s":mhd_src["s_mhd_Bz"]}
    components["psi-mhd"]   = {"f":mhd_var["f_mhd_psi"], "s":mhd_src["s_mhd_psi"]}

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

    #-------------------------------

    EXE="../../local/MFP.2d.gnu.MPI.ex"

    run = True


    res = [32, 64, 128] #, 256]
    max_step = 50

    data = {}

    L2_fail = 0
    L2_threshold = 0.95*2.0

    for r in res:
        name = "MMS%i"%r

        if run:
            inputs_name = make_inputs("MMS.inputs.template", name, {"nx":2*r, "ny":r, "name":name, "max_step":max_step})

            run = "rm -r {base_name}.chk*; mpirun -n 8 {exe} {name} 2>&1 | tee {name}.log".format(exe=EXE, base_name=name, name=inputs_name)
            print(run)

            subprocess.call(run, shell=True)

        files = get_files('.', include=[name, "chk%05i"%max_step], exclude=["input", ".png", ".py", ".lua", "old"], get_all=True)

        rh5 = ReadBoxLib(files[-1])

        # print("file = ",files[-1])

        t = rh5.time

        include = [] 
        for key, value in components.items():

            x, z = rh5.get("%s"%key)
            x, y = x
            
            if "sim_result" in value:
                value["sim_result"] += [z]
            else:
                value["sim_result"] = [z]

            y, x = np.meshgrid(y, x)

            if "x" in value:
                value["x"] += [x]
            else:
                value["x"] = [x]

            if "y" in value:
                value["y"] += [y]
            else:
                value["y"] = [y]
            
            f = value["f"](x, y, t)
            if type(f) != np.ndarray:
                f = np.ones(x.shape)*f
            
            if "mms_result" in value:
                value["mms_result"] += [f]
            else:
                value["mms_result"] = [f]

            s = value["s"](x, y, t)
            if type(s) != np.ndarray:
                s = np.ones(x.shape)*s

            if "mms_src" in value:
                value["mms_src"] += [s]
            else:
                value["mms_src"] = [s]

            err = z - f 
            L2 = np.linalg.norm(err,2)/float(x.size)

            if "err" in value:
                value["err"] += [err]
            else:
                value["err"] = [err]

            if "L2" in value:
                value["L2"] += [L2]
            else:
                value["L2"] = [L2]

            if "r" in value:
                value["r"] += [r]
            else:
                value["r"] = [r]

            if key not in include:
                    include.append(key)


    for key in include:
        data[key] = components[key]

    #-------------------------------



    axes = []

    nr = len(data)

    nc = int(np.sqrt(nr))
    nr = int(np.ceil(nr/float(nc)))

    fig = plt.figure(figsize=(6*nc,2*nr))

    gs = gridspec.GridSpec(ncols=3*nc, nrows=2*nr,
        width_ratios=nc*[1.0, 0.5, 0.5])

    items = [x for x in data.items()]
    N = len(items)
    cnt = 0
    for i in range(nr):
        for j in range(nc):

            if cnt > N - 1:
                continue
            
            key, dat = items[cnt]

            ax = fig.add_subplot(gs[2*i+0, 3*j+1]); axes.append(ax)
            ax.imshow(dat["sim_result"][-1].T, origin="lower")
            ax.set_title("sim")

            ax = fig.add_subplot(gs[2*i+1, 3*j+1]); axes.append(ax)
            ax.imshow(dat["mms_result"][-1].T, origin="lower")
            ax.set_title("exact")

            ax = fig.add_subplot(gs[2*i+0, 3*j+2]); axes.append(ax)
            ax.imshow(dat["err"][-1].T, origin="lower")
            ax.set_title("error")

            ax = fig.add_subplot(gs[2*i+1, 3*j+2]); axes.append(ax)
            ax.imshow(dat["mms_src"][-1].T, origin="lower")
            ax.set_title("src")

            ax = fig.add_subplot(gs[2*i:2*i+2, 3*j+0])
            l2 = dat["L2"]
            r = dat["r"]

            z = np.polyfit(np.log10(r),np.log10(l2),1)
            p = np.poly1d(z)

            if abs(z[0]) < L2_threshold:
                L2_fail = 1

            dat["O"] = z[0]

            ax.plot(r, 10**p(np.log10(r)), "r--", alpha=0.5, label="O(%0.3g)"%z[0])
            ax.loglog(r,l2,'kx')
            ax.legend()
            ax.set_title(key)

            cnt += 1

    for ax in axes:
        ax.set_aspect(1)
        ax.set_xticks([])
        ax.set_yticks([])    

    fig.tight_layout() #h_pad=0.05, w_pad=0.05)
    fig.savefig("plot.png",dpi=200)
    # plt.show()

    return L2_fail

if __name__ == "__main__":
    err = check()
    # print("MMS return code = ",err)
    sys.exit(err)
