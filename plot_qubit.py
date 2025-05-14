import matplotlib.pyplot as plt
import numpy as np

def read_file_columns(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            numbers = line.split()
            for i, number in enumerate(numbers):
                if len(data) <= i:
                    data.append([])
                data[i].append(float(number))
    return data

filein = open("params.txt")
tmax = float(filein.readline())

norms = read_file_columns("norms.txt")
OD = norms[0]
TD = norms[1]
incops = read_file_columns("sh_inc.txt")
sharpness = incops[0]
inc = incops[1]
inc_steer = incops[2]
funcs = read_file_columns("functions.txt")
lam = funcs[0]
b = funcs[1]

t = np.linspace(0, tmax, len(TD))

fig, ax = plt.subplots(1,2, figsize=(13,4), sharex=True, sharey=False, tight_layout=True)

ax[0].plot(t, OD, '-.', label=r'$D_\infty$', color='red')
ax[0].plot(t, TD, label=r'$D_1$', color='black')

size_inset = .35
axx = ax[0].inset_axes([0.05, 0.02, size_inset, size_inset])
axx.plot(t, lam, label=r'$\lambda$', color='r')
axx.plot(t, b, '--', label=r'$\beta$', color='b')
axx.set_xticks([])
axx.axvline(1, color='black', linewidth=1., linestyle='--')

ax[1].plot(t, inc, label=r'$I_{0,\text{steer}}$', color='black')
ax[1].plot(t, sharpness, '-.', label=r'$\Sigma$', color='red')

fontsize = 16

for i in range(2):
    ax[i].set_xlabel(r'$t$', fontsize=fontsize)
    ax[i].legend(loc = "upper left", fontsize=fontsize-2)
    ax[i].axvline(1, color='black', linewidth=1., linestyle='--')
    ax[i].tick_params(axis='both', which='major', labelsize=fontsize-3)
#ax[0].legend(loc = "upper right", fontsize=fontsize-2)
#ax[1].legend(loc = "upper left", fontsize=fontsize-2)
#ax[0].set_ylabel(r'$D_\infty(\Phi^*_t[X], \Phi^*_t[Y])$ and $D_1(\Phi_t[X], \Phi_t[Y])$', fontsize=fontsize)
ax[0].set_ylabel(r'$D_\infty$ and $D_1$', fontsize=fontsize)
ax[1].set_ylabel(r'$I_{0,\text{steer}}$ and $\Sigma$', fontsize=fontsize)
ax[0].set_ylim(-0.05, 1.05)
axx.legend(loc = "upper left", fontsize=fontsize-4)
axx.tick_params(axis='both', which='major', labelsize=fontsize-4)

plt.savefig('qubit.png', dpi=300)

plt.show()