import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.linalg import expm
from scipy.linalg import inv
from scipy.linalg import det
    
def a(tt):
    return 0
def b(tt):
    return 0
def c(tt):
    if tt < 1:
        return 0
    return tt - 1
    
def A(tt):
    return np.array([
        [0, a(tt), b(tt)],
        [-a(tt), 0, c(tt)],
        [-b(tt), -c(tt), 0]
    ])

def O(tt):
    return expm(A(tt))

lmin = .4
def l1 (tt):
    if tt < 1:
        return (1-lmin)*(1.-tt)*(1.-tt) + lmin
    else:
        return lmin   
def l2 (tt):
    return l1(tt) 
def l3 (tt):
    return 1

def D(tt):
    return np.array([
        [l1(tt), 0, 0],
        [0, l2(tt), 0],
        [0, 0, l3(tt)]
    ])

OD_or_DO = 'OD' # Set it to 'OD' or 'DO' depending on which evo want to compute
def Lambda(tt):
    if OD_or_DO == 'DO':
        return D(tt) @ inv(O(tt))
    return O(tt) @ D(tt)

# Generate initial unit sphere points
phi, theta = np.mgrid[0:np.pi:30j, 0:2*np.pi:60j]
x_sphere = np.sin(phi) * np.cos(theta)
y_sphere = np.sin(phi) * np.sin(theta)
z_sphere = np.cos(phi)

# Stack points into (3, N) shape for vectorized transformation
sphere_points = np.array([x_sphere.ravel(), y_sphere.ravel(), z_sphere.ravel()])

# Initialize figure
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

x_0, y_0, z_0 = sphere_points.reshape(3, *x_sphere.shape)

def set_axes_equal(ax):
    limits = np.array([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
    max_range = (limits[:, 1] - limits[:, 0]).max() / 2
    midpoints = np.mean(limits, axis=1)
    ax.set_xlim(midpoints[0] - max_range, midpoints[0] + max_range)
    ax.set_ylim(midpoints[1] - max_range, midpoints[1] + max_range)
    ax.set_zlim(midpoints[2] - max_range, midpoints[2] + max_range)

total_frames = 301
# Animation update function
def update(frame):
    t = frame / 100  # Time parameter
    transformed_points = Lambda(t) @ sphere_points  # Apply transformation
    
    # Reshape back into grid
    x_t, y_t, z_t = transformed_points.reshape(3, *x_sphere.shape)

    # Clear and update surface
    ax.clear()
    ax.plot_surface(x_t, y_t, z_t, color='c', alpha=0.9, edgecolor='black')
    ax.plot_surface(x_0, y_0, z_0, color='b', alpha=0.1)

    ax.plot([0, 1.4], [0, 0], [0, 0], 'black', linewidth=1)  # X-axis
    ax.plot([0, 0], [0, 1.4], [0, 0], 'black', linewidth=1)  # Y-axis
    ax.plot([0, 0], [0, 0], [0, 1.2], 'black', linewidth=1)  # Z-axis

    ax.text(1.6, 0, 0, r'$x$', color='black', fontsize=12)
    ax.text(0, 1.5, 0, r'$y$', color='black', fontsize=12)
    ax.text(0, 0, 1.3, r'$z$', color='black', fontsize=12)

    # Remove grid and box
    ax.set_axis_off()
    set_axes_equal(ax)

    ax.view_init(elev=15, azim=45)

    # Scale the z-axis to improve sphere proportions
    ax.set_box_aspect([1, 1, 1])  # Scaling factor for the z-axis

    if OD_or_DO == 'DO':
        ax.set_title(r'$D O^\top$, $t = $'+str(t))
    else:
        ax.set_title(r'$O D$, $t = $'+str(t))
    #if frame == total_frames - 1:
    #    plt.savefig(str(frame//100)+'_'+OD_or_DO+'.png', dpi=300)
    return ax,

# Create animation
ani = animation.FuncAnimation(fig, update, frames=total_frames, interval=50, blit=False)

# Save animation as video file
#if total_frames == 301:
#    ani.save(OD_or_DO+'.mp4', writer='ffmpeg', fps=30)

plt.show()