
import sys  # nopep8
cmd_folder = "../../../vis"  # nopep8
if cmd_folder not in sys.path:  # nopep8
    sys.path.insert(0, cmd_folder)

from tile_mov import tile_movie
from make_mov import make_all
import numpy as np


# ==============================================================================
# MAKE MOVIES
# ==============================================================================

if 1:
    dt = 0.005

    Q = []

    q = {}
    q["h5_dir"] = "."
    q["level"] = -1

    q["select"] = []
    q["label"] = []
    q["cmap"] = []
    q["scale_limits"] = []
    q["label_colour"] = []

    # HYDRO
    q["select"] += [{"exp": "{density-ion}", "tag": "density-ion", "label": r"$\rho_i$",
                     "cmap": {"map": "magma_r", "centered": False}, "scale_limits": 1.0, "label_colour": 'white',
                     "particles":[{"select":"ion", "style":{"s":0.5, "c":"k", "alpha":0.3}}],
                     "grid":{"color":"w", "opacity":1.0, "linewidth":0.1}
                     }]

    q["select"] += [{"exp": "{density-electron}", "tag": "density-electron", "label": r"$\rho_e$",
                     "cmap": {"map": "magma_r", "centered": False}, "scale_limits": 1.0, "label_colour": 'white',
                     "particles":[{"select":"electron", "streak":{"length":25, "style":{"c":"w", "lw":0.5, "alpha":0.3}}, "style":{"s":0.5, "c":"k", "alpha":0.3}
                     }]}]

    ##
    q["framerate"] = 20
    q["mov_save"] = q["h5_dir"] + "/mov"
    q["minmax"] = False
    q["time"] = True
    q["colorbar"] = True
    q["insert_label"] = True
    q["text_size"] = 12
    q["offset"] = [0.0, 0.0]
    q["xy_limits"] = [[-0.5, 1.5], [0.0, 1.0]]
    # q["xy_limits_final"] =
    q["file_include"] = ["TRMI.plt"]
    q["file_exclude"] = ["chk"]
    q["cores"] = 10
    q["frame_size"] = (6, 3)
    q["time_span"] = [] #np.arange(0,0.1,0.005).tolist()
    q["force_data"] = False
    q["force_frames"] = True
    q["only_frames"] = False
    q["dpi"] = 200

    q["normalize"] = "all"

    Q.append(q)

    make_all(Q)

# ==============================================================================
# TILE MOVIES
# ==============================================================================

if 1:
    q = {}
    q["dirs"] = []
    q["dirs"].append({"dir": "./mov/density-ion",
                      "level": -1, "tile": [0, 0], "label": ""})
    q["dirs"].append({"dir": "./mov/density-electron",
                      "level": -1, "tile": [1, 0], "label": ""})

    q["name"] = "TRMI"

    q["out_dir"] = "./mov"

    q["trim"] = [0, 999999]
    q["crop"] = False  # [[180,1150],[0,530]]
    q["max_individual_size"] = (1920, 1080)
    q["label_size"] = 16
    q["framerate"] = 30
    q["cores"] = 4
    q["level_label"] = False

    tile_movie(q)

print("DONE")
