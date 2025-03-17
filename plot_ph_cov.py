import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Reading the rates
filein = open("rates.dat")
tf = float(filein.readline())
dt = float(filein.readline())

Npoints = int(tf/dt)
t_R = np.arange(0,tf,dt)

gamma_p = np.zeros(Npoints)
gamma_m = np.zeros(Npoints)
gamma_z = np.zeros(Npoints)
xi_p = np.zeros(Npoints)
xi_m = np.zeros(Npoints)
xi_z = np.zeros(Npoints)

for i in range(Npoints):
    rates = filein.readline().split()
    gamma_p[i] = rates[0]
    gamma_m[i] = rates[1]
    gamma_z[i] = rates[2]
    rates = filein.readline().split()
    xi_p[i] = rates[0]
    xi_m[i] = rates[1]
    xi_z[i] = rates[2]

# Reading the norms
filein = open("norm.dat")
tf = float(filein.readline())
dt = float(filein.readline())

Npoints = int(tf/dt)
t_N = np.arange(0,tf,dt)

E = np.zeros(Npoints)
F = np.zeros(Npoints)
E_F = np.zeros(Npoints)

for i in range(Npoints):
    norms = filein.readline().split()
    E[i] = norms[0]
    F[i] = norms[1]
    E_F[i] = norms[2]

fig, ax = plt.subplots(1,2, figsize=(13,5), sharex=False, sharey=False, tight_layout=False)

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
ax[0].plot(t_R, gamma_p, label = r'$\gamma_+$', color=colors[0])
ax[0].plot(t_R, gamma_m, label = r'$\gamma_-$', color=colors[1])
ax[0].plot(t_R, gamma_z, label = r'$\gamma_z$', color=colors[2])
ax[0].plot(t_R, xi_p, '--', label = r'$\xi_+$', color=colors[0])
ax[0].plot(t_R, xi_m, '--', label = r'$\xi_-$', color=colors[1])
ax[0].plot(t_R, xi_z, '--', label = r'$\xi_z$', color=colors[2])
ax[0].axhline(0, color='black', linewidth=1.)

ax[1].plot(t_N, E, label = r'$E$')
ax[1].plot(t_N, F, label = r'$F$')
ax[1].plot(t_N, E_F, label = r'$E-F$')

for i in range(2):
    ax[i].legend(loc = "lower left")
    ax[i].set_xlabel(r'$t$')
ax[0].set_title('Rates')
ax[1].set_title(r'$\Vert \Phi^\ddag_t[X] \Vert_\infty$')

plt.savefig('ph_cov.png', dpi=300)

plt.show()