
import sympy as sp
from sympy import Rational, Integer, symbols, simplify, latex, sqrt, pi, diff, sin, cos, exp
from sympy import init_printing, pprint
from sympy.utilities.lambdify import lambdify
from numpy.random import random

init_printing(num_columns=1000)

zero = Integer(0)
half = Rational(1,2)
one = Integer(1)
two = Integer(2)

#==============================================================================
# MMS sources
#==============================================================================

x, y, t = symbols('x y t', real=True)

gamma = Rational(5,3)
c0 = Integer(10)
dL = half
dD = half

cr = 0.5

ch_maxwell = c0
ch_mhd = 1.5

c_damp_maxwell = ch_maxwell/cr
c_damp_mhd = ch_mhd/cr


hydro_const = [{"q":one, "m":one, "g":gamma}, {"q":-one, "m":0.1, "g":gamma}]
field_const = {"c0":c0, "dL":dL, "dD":dD, "ch_maxwell":ch_maxwell, "cr":cr}
mhd_const = {"mhd_g":gamma, "mhd_m":one, "ch_mhd":ch_mhd}

def MMS(x, y, t, a0, ax, ay, axx=2*pi, ayy=2*pi, L=pi):
    return a0 + ax*sin(axx*pi*x/L + 10*t) + ay*cos(ayy*pi*y/L + 5*t)
    # return a0 + ax*sin(axx*pi*x/L) + ay*cos(ayy*pi*y/L)

rho = MMS(x,y,t, 1.0, 0.44, 0.10)
u   = MMS(x,y,t, 0.0, 0.51, 0.23)
v   = MMS(x,y,t, 0.0, 0.60, 0.81)
w   = MMS(x,y,t, 0.0, 0.27, 0.13)
p   = MMS(x,y,t, 1.0, 0.20, 0.30)

Bx  = MMS(x,y,t, 0.0, 0.93, 0.39)
By  = MMS(x,y,t, 0.0, 0.17, 0.89)
Bz  = MMS(x,y,t, 0.0, 0.98, 0.38)
psi = MMS(x,y,t, 0.0, 0.23, 0.53)

Dx  = MMS(x,y,t, 0.0, 0.28, 0.81)
Dy  = MMS(x,y,t, 0.0, 0.57, 0.10)
Dz  = MMS(x,y,t, 0.0, 0.32, 0.70)
phi = MMS(x,y,t, 0.0, 0.21, 0.35)

Ex = Dx
Ey = Dy
Ez = Dz

# hydro species

hydro_var = []
hydro_src = []

jx = 0
jy = 0
jz = 0
cd = 0

for n in range(len(hydro_const)):

    q = hydro_const[n]["q"]
    m = hydro_const[n]["m"]
    g = hydro_const[n]["g"]

    mx = u*rho
    my = v*rho
    mz = w*rho
    nrg = p/(g - one) + half*rho*(u**2 + v**2 + w**2)

    H = {
        "rho":rho, "u":u, "v":v, "w":w, "p":p,
        "mx":mx, "my":my, "mz":mz, "nrg":nrg,
    }

    S = {}

    # continuity
    drho_dt = diff(rho, t)

    div_rho_dx = diff(mx, x)
    div_rho_dy = diff(my, y)

    S["rho"] = drho_dt + div_rho_dx + div_rho_dy

    cd += q*rho/m

    # x momentum

    dmx_dt = diff(mx, t)

    div_mx_x = diff(rho*u*u + p, x)
    div_mx_y = diff(rho*u*v, y)

    s_mx = (rho*q)/(m*dL)*(c0*Ex + v*Bz - w*By)

    S["mx"] = dmx_dt + div_mx_x + div_mx_y - s_mx

    jx += u*q*rho/m

    # y momentum

    dmy_dt = diff(my, t)

    div_my_x = diff(rho*u*v, x)
    div_my_y = diff(rho*v*v + p, y)

    s_my = (rho*q)/(m*dL)*(c0*Ey + w*Bx - u*Bz)

    S["my"] = dmy_dt + div_my_x + div_my_y - s_my

    jy += v*q*rho/m

    # z momentum

    dmz_dt = diff(mz, t)

    div_mz_x = diff(rho*u*w, x)
    div_mz_y = diff(rho*v*w, y)

    s_mz = (rho*q)/(m*dL)*(c0*Ez + u*By - v*Bx)

    S["mz"] = dmz_dt + div_mz_x + div_mz_y - s_mz

    jz += w*q*rho/m

    # energy

    dnrg_dt = diff(nrg, t)

    div_nrg_x = diff(u*(nrg + p), x)
    div_nrg_y = diff(v*(nrg + p), y)

    s_nrg = (q*c0)/(m*dL)*(mx*Ex + my*Ey + mz*Ez)

    S["nrg"] = dnrg_dt + div_nrg_x + div_nrg_y - s_nrg

    hydro_var.append(H)
    hydro_src.append(S)

field_var = {"Bx":Bx, "By":By, "Bz":Bz, "Dx":Dx, "Dy":Dy, "Dz":Dz, "phi":phi, "psi":psi}
field_src = {}

d_Dx_t = diff(Dx, t)
div_Dx_x = c0*diff(phi,x)
div_Dx_y = -c0*diff(Bz,y)
field_src["Dx"] = d_Dx_t + div_Dx_x + div_Dx_y + dL/(c0*dD**2)*jx

d_Dy_t = diff(Dy, t)
div_Dy_x = c0*diff(Bz,x)
div_Dy_y = c0*diff(phi,y)
field_src["Dy"] = d_Dy_t + div_Dy_x + div_Dy_y + dL/(c0*dD**2)*jy

d_Dz_t = diff(Dz, t)
div_Dz_x = -c0*diff(By,x)
div_Dz_y =  c0*diff(Bx,y)
field_src["Dz"] = d_Dz_t + div_Dz_x + div_Dz_y + dL/(c0*dD**2)*jz

d_Bx_t = diff(Bx,t)
div_Bx_x = c0*diff(psi,x)
div_Bx_y = c0*diff(Dz,y)
field_src["Bx"] = d_Bx_t + div_Bx_x + div_Bx_y

d_By_t = diff(By,t)
div_By_x = -c0*diff(Dz,x)
div_By_y = c0*diff(psi,y)
field_src["By"] = d_By_t + div_By_x + div_By_y

d_Bz_t = diff(Bz,t)
div_Bz_x = c0*diff(Dy,x)
div_Bz_y = -c0*diff(Dx,y)
field_src["Bz"] = d_Bz_t + div_Bz_x + div_Bz_y

d_phi_t = diff(phi,t)
div_phi_x = (ch_maxwell**2/c0)*diff(Dx,x)
div_phi_y = (ch_maxwell**2/c0)*diff(Dy,y)
field_src["phi"] = d_phi_t + div_phi_x + div_phi_y - (dL*ch_maxwell**2)/(c0**2*dD**2)*cd + c_damp_maxwell*phi

d_psi_t = diff(psi,t)
div_psi_x = (ch_maxwell**2/c0)*diff(Bx,x)
div_psi_y = (ch_maxwell**2/c0)*diff(By,y)
field_src["psi"] = d_psi_t + div_psi_x + div_psi_y + c_damp_maxwell*psi

## MHD

nrg = p/(g - one) + half*(rho*(u**2 + v**2 + w**2) + Bx**2 + By**2 + Bz**2)

mhd_var = {
    "mhd_rho":rho, "mhd_u":u, "mhd_v":v, "mhd_w":w, "mhd_p":p,
    "mhd_mx":mx, "mhd_my":my, "mhd_mz":mz, "mhd_nrg":nrg,
    "mhd_Bx":Bx, "mhd_By":By, "mhd_Bz":Bz, "mhd_psi":psi
}

mhd_src = {}

# continuity
drho_dt = diff(rho, t)
div_rho_dx = diff(mx, x)
div_rho_dy = diff(my, y)

mhd_src["mhd_rho"] = drho_dt + div_rho_dx + div_rho_dy


# x momentum

dmx_dt = diff(mx, t)

div_mx_x = diff(rho*u*u + p + half*(By*By + Bz*Bz - Bx*Bx), x)
div_mx_y = diff(rho*u*v - Bx*By, y)

mhd_src["mhd_mx"] = dmx_dt + div_mx_x + div_mx_y

# y momentum

dmy_dt = diff(my, t)

div_my_x = diff(rho*u*v - Bx*By, x)
div_my_y = diff(rho*v*v + p + half*(Bx*Bx + Bz*Bz - By*By), y)

mhd_src["mhd_my"] = dmy_dt + div_my_x + div_my_y

# z momentum

dmz_dt = diff(mz, t)

div_mz_x = diff(rho*u*w - Bx*Bz, x)
div_mz_y = diff(rho*v*w - By*Bz, y)

mhd_src["mhd_mz"] = dmz_dt + div_mz_x + div_mz_y

# energy

dnrg_dt = diff(nrg, t)

div_nrg_x = diff((nrg + p + half*(Bx*Bx + By*By + Bz*Bz))*u - Bx*(u*Bx + v*By + w*Bz),x)
div_nrg_y = diff((nrg + p + half*(Bx*Bx + By*By + Bz*Bz))*v - By*(u*Bx + v*By + w*Bz),y)

mhd_src["mhd_nrg"] = dnrg_dt + div_nrg_x + div_nrg_y

# Bx

dBx_dt = diff(Bx, t)

div_Bx_x = diff(psi,x)
div_Bx_y = diff(v*Bx - u*By,y)

mhd_src["mhd_Bx"] = dBx_dt + div_Bx_x + div_Bx_y

# By

dBy_dt = diff(By, t)

div_By_x = diff(u*By - v*Bx,x)
div_By_y = diff(psi,y)

mhd_src["mhd_By"] = dBy_dt + div_By_x + div_By_y

# Bz

dBz_dt = diff(Bz, t)

div_Bz_x = diff(u*Bz - w*Bx,x)
div_Bz_y = diff(v*Bz - w*By,y)

mhd_src["mhd_Bz"] = dBz_dt + div_Bz_x + div_Bz_y

# psi

dpsi_dt = diff(psi, t)

div_psi_x = diff(ch_mhd**2*Bx,x)
div_psi_y = diff(ch_mhd**2*By,y)

mhd_src["mhd_psi"] = dpsi_dt + div_psi_x + div_psi_y + c_damp_mhd*psi

#==============================================================================
# output
#==============================================================================

def function_template(name, guts):

    for v in ["pi", "sin", "cos", "exp"]:
        guts = guts.replace(v, "math.%s"%v)

    guts = guts.replace("**", "^")

    s =  "function %s(dat)\n"%name
    s += "  x = dat['x']\n"
    s += "  y = dat['y']\n"
    s += "  z = dat['z']\n"
    s += "  t = dat['t']\n"
    s += "  return %s\n"%guts
    s += "end\n"

    return s

fid = open("mms.lua", "w")

for i, cst in enumerate(hydro_const):
    for key, val in cst.items():
        fid.write("mms_%s_%i = %s\n"%(key,i,str(val)))

for key, val in field_const.items():
    fid.write("mms_%s = %s\n"%(key,str(val)))

for key, val in mhd_const.items():
    fid.write("mms_%s = %s\n"%(key,str(val)))

for i, var in enumerate(hydro_var):
    func = {}
    for key, val in var.items():
        func["f_%s_%i"%(key,i)] = lambdify((x, y, t), val)
        fid.write(function_template("mms_%s_%i"%(key,i), str(val)))
    var.update(func)

for i, src in enumerate(hydro_src):
    func = {}
    for key, val in src.items():
        func["s_%s_%i"%(key,i)] = lambdify((x, y, t), val)
        fid.write(function_template("mms_S%s_%i"%(key,i), str(val)))
    src.update(func)


func = {}
for key, val in field_var.items():
    func["f_%s"%(key)] = lambdify((x, y, t), val)
    fid.write(function_template("mms_%s"%(key), str(val)))
field_var.update(func)

func = {}
for key, val in field_src.items():
    func["s_%s"%(key)] = lambdify((x, y, t), val)
    fid.write(function_template("mms_S%s"%(key), str(val)))
field_src.update(func)

func = {}
for key, val in mhd_var.items():
    func["f_%s"%(key)] = lambdify((x, y, t), val)
    fid.write(function_template("mms_%s"%(key), str(val)))
mhd_var.update(func)

func = {}
for key, val in mhd_src.items():
    func["s_%s"%(key)] = lambdify((x, y, t), val)
    fid.write(function_template("mms_S%s"%(key), str(val)))
mhd_src.update(func)

fid.close()

