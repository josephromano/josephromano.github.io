import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.special as special

phi0 = 40 # max angular displacement (degrees)

# convert to radians
phi0 = np.deg2rad(phi0)

# some constants
g = 9.8 # m/s^2
l = 1. # m
omega0 = np.sqrt(g/l)
P0 = 2*np.pi/omega0
	
# eccentricity parameter for elliptic functions
k = np.sin(phi0/2)
m = k**2

# discrete times	
t = np.linspace(0, 4*P0, 10000)

# small angle approx
phi_approx = phi0 *np.cos(omega0*t)
	
# exact solution
x = omega0*t
P = (4/omega0)*special.ellipk(m)
sn_x, cn_x, dn_x, ph_x = special.ellipj(x+omega0*P/4, m)
phi = 2*np.arcsin(k*sn_x)

# define figure
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-1.25, 1.25), ylim=(-1.25, 1.25))
ax.grid()
line1, = ax.plot([], [], 'o--', color='grey', lw=2)
line2, = ax.plot([], [], 'o-', color='black', lw=2)
ax.text(-1,0.5, 'grey, dashed: small-angle approx')
ax.text(-1,0.75, 'black, solid: exact solution')

def update(num, phi_approx, phi, line1, line2):
	line1.set_data([0,np.sin(phi_approx[num])],[0,-np.cos(phi_approx[num])])
	line2.set_data([0, np.sin(phi[num])], [0,-np.cos(phi[num])])
	return line1, line2

#def init():
#	line.set_data([0, np.sin(phi_approx[0])], [0,-np.cos(phi_approx[0])])
#	return line

ani = animation.FuncAnimation(fig, update, len(t), fargs=[phi_approx,
phi, line1, line2], interval=1, blit=True, repeat=False)

plt.show()

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# writer = animation.FFMpegWriter(
#     fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)	
