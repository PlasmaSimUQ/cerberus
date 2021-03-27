import numpy as np
from scipy.special import jv, hankel2
import pylab as plt

def Ez(x,y,a=0.5,k=2*np.pi, rho_d=2.0, phi_d=0.0, S=20):

    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)

    E0 = 0
    E1 = 0

    for n in range(-S,S+1):

        c = - jv(n, k*a)/hankel2(n, k*a)

        E0 += hankel2(n, k*rho_d)*(jv(n, k*rho) + c*hankel2(n, k*rho))*np.exp(1j*n*(phi-phi_d))
        E1 += hankel2(n, k*rho)*(jv(n, k*rho_d) + c*hankel2(n, k*rho_d))*np.exp(1j*n*(phi-phi_d))

    mask = rho <= rho_d

    E1[mask] = E0[mask]

    return E1

N = 200
x = np.linspace(-4,6,N)
y = np.linspace(-4,4,N)

x, y = np.meshgrid(x,y)

r = np.sqrt(x**2 + y**2)

mask = r <= 0.5

x = np.ma.masked_where(mask,x)
y = np.ma.masked_where(mask,y)

ez = Ez(x,y)

fig = plt.figure()

ax = fig.add_subplot(1,1,1)
ax.contourf(x,y,ez)
ax.set_aspect(1)

plt.show()
