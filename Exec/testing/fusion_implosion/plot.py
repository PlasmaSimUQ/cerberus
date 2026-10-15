import os
import re
import sys

cmd_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../vis")
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files


class ReadBoxLibHe4(ReadBoxLib):
    """ReadBoxLib with a derived 'He4' component: the ash mass fraction is the
    implicit last ion component, 1 - alpha_D - alpha_T."""

    def get(self, name, **kw):
        if name != "He4":
            return super().get(name, **kw)
        xs, aD = super().get("alpha_0-ion", **kw)
        _, aT = super().get("alpha_1-ion", **kw)
        return xs, 1.0 - aD - aT


ReadBoxLib = ReadBoxLibHe4

import numpy as np

# vis/get_boxlib.py still uses np.product, removed in numpy 2.0
if not hasattr(np, "product"):
    np.product = np.prod

import matplotlib

matplotlib.use("Agg")
from matplotlib.colors import LogNorm, Normalize
import pylab as plt
from matplotlib.image import NonUniformImage

# plotfile prefix: default for the top-level case, or the case name when run
# from one of the comparison subdirectories (python3 ../plot.py <case>)
plt_file = sys.argv[1] if len(sys.argv) > 1 else "fusion_implosion"
fig_dir = "figures"

# the fields drawn, one panel each: (plotfile component, label, cmap, log?)
FIELDS = [
    ("rho-ion", r"$\rho_i$", "viridis", True),
    ("T-ion", r"$T_i$ (code)", "magma", True),
    ("T-electron", r"$T_e$ (code)", "magma", True),
    ("He4", r"He4 mass fraction", "coolwarm", False),
]

files = get_files(
    '.',
    include=[plt_file],
    # .chk holds only conserved variables, so it must not be picked up here
    exclude=["temp", ".chk", ".png", ".inputs", ".txt", ".py"],
    times=[],
    tol=1e-4,
    get_all=True,
)

if not files:
    raise SystemExit("no %s.plt* files found -- run the case first" % plt_file)

os.makedirs(fig_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# pass 1: colour limits, held FIXED across every frame so that the sequence is
# comparable frame-to-frame (a per-frame autoscale hides the compression on the
# way in and the decay on the way out)
# ---------------------------------------------------------------------------
clim = {name: [np.inf, -np.inf] for name, _, _, _ in FIELDS}

for fname in files:
    rh5 = ReadBoxLib(fname, max_level=-1)
    for name, _, _, _ in FIELDS:
        _, data = rh5.get(name)
        clim[name][0] = min(clim[name][0], np.min(data))
        clim[name][1] = max(clim[name][1], np.max(data))

# the ash fraction starts at exactly 0 (round-off negative is possible); pin
# its floor so the colour scale starts at zero
clim["He4"][0] = 0.0

# ---------------------------------------------------------------------------
# pass 2: one png per plotfile
# ---------------------------------------------------------------------------
def frame_tag(fname):
    """the step number encoded in RMI_implosion.plt00123"""
    match = re.search(r"\.plt(\d+)", os.path.basename(fname))
    return match.group(1) if match else os.path.basename(fname)


def draw(fname, out_name, title):
    rh5 = ReadBoxLib(fname, max_level=-1)

    fig, axes = plt.subplots(1, len(FIELDS), figsize=(5 * len(FIELDS), 5))

    for ax, (name, label, cmap, log) in zip(np.atleast_1d(axes), FIELDS):
        (xs, ys), data = rh5.get(name)

        window = [[xs[0], xs[-1]], [ys[0], ys[-1]]]
        vmin, vmax = clim[name]
        # density and pressure both span more than an order of magnitude across
        # the four layers -> log scale keeps the target/shell contrast visible
        # alongside the thin background.  NonUniformImage takes a norm, not
        # vmin/vmax, so the linear case gets an explicit Normalize.
        norm = LogNorm(vmin=vmin, vmax=vmax) if log else Normalize(vmin=vmin, vmax=vmax)

        im = NonUniformImage(
            ax, interpolation='nearest', extent=np.ravel(window), cmap=cmap, norm=norm
        )
        im.set_data(xs, ys, data.T)
        ax.add_image(im)
        ax.set_xlim(window[0])
        ax.set_ylim(window[1])
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.set_aspect(1)
        plt.colorbar(im, ax=ax, label=label, orientation="horizontal")

    fig.suptitle("%s  (t = %g)" % (title, rh5.time))
    fig.tight_layout()
    fig.savefig(out_name, dpi=200)
    plt.close(fig)


for fname in files:
    tag = frame_tag(fname)
    out_name = os.path.join(fig_dir, "%s_%s.png" % (plt_file, tag))
    draw(fname, out_name, "step %s" % tag.lstrip("0").rjust(1, "0"))
    print("wrote " + out_name)

# first and last available plotfile at the top level too, so the initial
# condition can be checked without digging through figures/
for tag, fname in [("initial", files[0]), ("final", files[-1])]:
    out_name = "%s_%s.png" % (plt_file, tag)
    draw(fname, out_name, tag)
    print("wrote " + out_name)

print("DONE")
