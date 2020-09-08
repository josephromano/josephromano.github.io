import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.optimize as optimize

# numerically solve for the motion of two masses connected by a 
# string passing through a hole on a frictionless table, in a 
# uniform gravitational field g

# set parameters
g = 9.8 # m/s^2
#m1 = 0.5 # kg (on table) (closed orbit)
m1 = 0.65 # kg (on table)
m2 = 1.0 # kg (below table)
r0 = 0.5 # m  (radius for stable circular orbit)
Mz = np.sqrt(m1*m2*g*r0**3) # (angular momentum needed for r=r0)
omega_phi = np.sqrt((m2/m1)*(g/r0)) # angular frequency for a circular orbit
P_phi = 2*np.pi/omega_phi # corresponding angular period
 
# discrete times, radii, phi
N = 4*1000
t = np.linspace(0, 4*P_phi, N)
dt = t[1]-t[0]
r = np.zeros(N)
phi = np.zeros(N)

# initial value of r, phi
rmin = 0.7*r0 # minimum r value
r[0] = rmin
phi[0] = 0

# total energy for the system
E = m2*g*rmin + Mz**2/(2*m1*rmin**2)
print('E =', E)

# determine rmax by solving cubic equation
def f(rmax, m1, m2, g, Mz):
	y = E - ( m2*g*rmax + Mz**2/(2*m1*rmax**2) )
	return y

rmax = optimize.newton(f, 1.1*r0, args=(m1, m2, g, Mz))
print('rmin =', rmin)
print('r0 =', r0)
print('rmax =', rmax)

# plot effective potential
rs = np.linspace(0.2, 1.2*rmax, 1000)
Ueff = m2*g*rs + Mz**2/(2*m1*rs**2)
plt.figure()
plt.plot(rs, Ueff)
plt.axhline(E, color='grey', ls='dashed')
plt.axvline(rmin, color='grey')
plt.axvline(rmax, color='grey')
plt.xlabel('r')
plt.ylabel('Ueff')

# integrate equations of motion
sign = 1
for ii in range(1,len(t)):

    # radial equation
	dr = sign * dt*np.sqrt((2/(m1+m2))* (E- m2*g*r[ii-1]-Mz**2/(2*m1*r[ii-1]**2)))
	r[ii] = r[ii-1]+dr
	if r[ii] > rmax:
		r[ii] = r[ii]-abs(dr)
		sign = -sign
	if r[ii] < rmin:
		r[ii] = r[ii]+abs(dr)
		sign = -sign
	
	# angular equation (from conservation of angular momentum)
	dphi = dt*Mz/(m1*r[ii-1]**2)
	phi[ii] = phi[ii-1]+dphi

# plot r vs. phi
plt.figure()
plt.plot(phi, r)
plt.xlabel('phi')
plt.ylabel('r')

# cartesian coordinates
x = r*np.cos(phi)
y = r*np.sin(phi)

# plot x vs. y
plt.figure()
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-1,1)
plt.ylim(-1,1)

# define figure
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-1., 1.), ylim=(-1., 1.))
ax.grid()
line, = ax.plot([], [], 'o-', color='grey', lw=2)

def update(num, x, y, line):
	line.set_data([0,x[num]],[0,y[num]])
	return line,

ani = animation.FuncAnimation(fig, update, len(phi), fargs=[x, y, line], interval=1, blit=True, repeat=False)

plt.show()

# old code that doesn't using animation function
#for ii in range(N):     
#	line.set_data([0, x[ii]],[0,y[ii]])
#	plt.pause(0.001)

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# writer = animation.FFMpegWriter(
#     fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)
##############################################
