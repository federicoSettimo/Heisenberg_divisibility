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

fig, ax = plt.subplots(1,1, figsize=(9,5), sharex=False, sharey=False, tight_layout=False)

ax.plot(t_N, E, '--', label = r'$X=\vert0\rangle\langle0\vert$, $Y=0$', color='green')
ax.plot(t_N, F, '-.', label = r'$X=0$, $Y=\vert1\rangle\langle1\vert$', color='red')
ax.plot(t_N, E_F, label = r'$X=\vert0\rangle\langle0\vert$, $Y=\vert1\rangle\langle1\vert$', color='black')

colors = ['r', 'b']
axx = inset_axes(ax, width="33%", height="33%", loc=1)
axx.plot(t_R, gamma_p, label = r'$\gamma_\pm$', color=colors[0])
axx.plot(t_R, gamma_m, '--', color=colors[0])
axx.plot(t_R, xi_p, label = r'$\xi_\pm$', color=colors[1])
axx.plot(t_R, xi_m, '--', color=colors[1])
axx.axhline(0, color='black', linewidth=1.)

t_NM = 1.88
ax.axvline(t_NM, color='black', linewidth=1., linestyle='--')
axx.axvline(t_NM, color='black', linewidth=1., linestyle='--')


fontsize = 16
ax.legend(loc = "lower left", fontsize=fontsize-2)
ax.set_xlabel(r'$t$', fontsize=fontsize)
ax.set_ylabel(r'$D_\infty(\Phi^*_t[X], \Phi^*_t[Y])$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=fontsize-3)

axx.legend(loc = "lower left", fontsize=fontsize-4)
axx.set_xlabel(r'$t$', fontsize=fontsize-2)
axx.tick_params(axis='both', which='major', labelsize=fontsize-5)

plt.savefig('ph_cov.png', dpi=300)

plt.show()