import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # noqa F401


def kwargs_plot_cylinder( r, h, n=1000 ):
    z = np.asarray([0.0,h],dtype=float)
    theta = np.linspace(0, 2 * np.pi, n)
    theta_grid, z_grid = np.meshgrid(theta, z)
    #x_grid = r * np.cos(theta_grid)
    #y_grid = r * np.sin(theta_grid)
    return dict( X = r * np.cos(theta_grid),
                 Y = r * np.sin(theta_grid),
                 Z = z_grid )

def kwargs_plot_sphere( r, n=1000 ):
    u = np.linspace(0, 2 * np.pi, n)  # azimuthal angle
    v = np.linspace(0, np.pi, n)       # polar angle
    return dict( X = r * np.outer(np.cos(u), np.sin(v)),
                 Y = r * np.outer(np.sin(u), np.sin(v)),
                 Z = r * np.outer(np.ones(np.size(u)), np.cos(v)) )

# Cylinder parameters
r = 1  # radius
h = 3  # height
n_points = 1000  # Resolution
#
## Create the cylinder
#z_cylinder = np.linspace(0, h, 2)#n_points)
#theta = np.linspace(0, 2 * np.pi, n_points)
#theta_grid, z_grid = np.meshgrid(theta, z_cylinder)
#x_grid = r * np.cos(theta_grid)
#y_grid = r * np.sin(theta_grid)
#
# Create the figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the cylinder surface
ax.plot_surface( **kwargs_plot_sphere(r=r,n=1000),
                 alpha=0.5, color='green' )
#ax.plot_surface( **kwargs_plot_cylinder(r=r,h=h,n=1000),
#                 alpha=0.5, color='green' )
#ax.plot_surface(x_grid, y_grid, z_grid, alpha=0.5, color='cyan')

# Define rays intersecting the cylinder (random for illustration)
num_rays = 4
ray_angles = np.linspace(0, 2 * np.pi, num_rays, endpoint=False)
for angle in ray_angles:
    # Ray from above hitting the cylinder
    x_ray = np.array([r * np.cos(angle), r * np.cos(angle)])
    y_ray = np.array([r * np.sin(angle), r * np.sin(angle)])
    z_ray = np.array([h, -1])  # Starts above the cylinder

    # Plot the rays
    ax.plot(x_ray, y_ray, z_ray, color='red', linewidth=2)

# Set labels and view angle
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.view_init(elev=30, azim=30)

plt.title('Cylinder with Rays Hitting It')
plt.show()
