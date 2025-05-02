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

fig, ax = plt.subplots(1,2, figsize=(13,4), sharex=True, sharey=True, tight_layout=True)

ax[0].plot(t, OD, label=r'$D_\infty$', color='red')
ax[0].plot(t, TD, '--', label=r'$D_1$', color='blue')

size_inset = .3
axx = ax[0].inset_axes([0.98-size_inset, 0.02, size_inset, size_inset])
axx.plot(t, lam, label=r'$\lambda$')
axx.plot(t, b, label=r'$\beta$')
axx.set_xticks([])

ax[1].plot(t, inc, label=r'$I_0$', color='red')
ax[1].plot(t, inc_steer, label=r'$I_{\text{steer}}$', color='blue')
ax[1].plot(t, sharpness, '--', label=r'$\Sigma$', color='darkgreen')
#axx_t = ax[1].twinx()
#axx_t.plot(t, sharpness, '--', label=r'$\Sigma$', color='darkgreen')

for i in range(2):
    ax[i].set_xlabel(r'$t$')
    ax[i].legend(loc = "lower left")
    ax[i].axvline(1, color='red', linewidth=.5)
ax[0].set_title(r'$D_\infty(\Phi^*_t[X], \Phi^*_t[Y])$ and $D_1(\Phi_t[X], \Phi_t[Y])$')
ax[1].set_title(r'Incompatibility and sharpness')
axx.legend(loc = "lower right")

plt.show()