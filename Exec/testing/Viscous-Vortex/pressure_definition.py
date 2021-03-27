import numpy as np
from scipy.special import expi

G = 1.0
t0 = 10.0
mu0 = 0.08
rho0 = 1.0
p0 = 1.00548676
nu = mu0/rho0

def P(r):
    return p0 - (G**2*rho0*((-1 + np.exp(r**2/(4.*t0*nu)))**2/(2.*np.exp(r**2/(2.*t0*nu))*r**2) + \
    (expi(-r**2/(2.*t0*nu)) - expi(-r**2/(4.*t0*nu)))/(4.*t0*nu)))/(4.*np.pi**2)

def make_lua_data(name, x, y):

    f = open(name + ".lua", "w")

    f.write("x_data_%s = {"%name)
    for i, x_ in enumerate(x):
        f.write(str(x_))
        if i < len(x) - 1:
            f.write(",")
        if not i%5:
            f.write("\n")
    f.write("}\n")

    f.write("\n")

    f.write("y_data_%s = {"%name)
    for i, y_ in enumerate(y):
        f.write(str(y_))
        if i < len(y) - 1:
            f.write(",")
        if not i%5:
            f.write("\n")
    f.write("}\n")

    f.close()


if __name__ == "__main__":

    N = 1000
    r = np.linspace(1e-6, 30, N)
    p = P(r)

    make_lua_data("pressure", r, p)

    if 0:
        import pylab as plt
        plt.plot(r,p)
        plt.show()
