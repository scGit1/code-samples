#!/usr/bin/python
#
# script to compute 2d wave propagation using the leap frog technique (numerical integration)

import numpy as np
import matplotlib.pyplot as plt
cmap = plt.cm.get_cmap('seismic')
cmap.set_under("yellow")
cmap.set_over("yellow")
#
# set up grid...including dt
# note that when setting 
# up the meshgrid, y, x = mgrid
# this makes the indexing of the 
# field be u[x_i,y_i]
#
N   = 150 ; L  = 1.
ds  = 2*L/float(N)
xy  = np.linspace(-L,L,N+2)
y,x = np.meshgrid(xy,xy)
dt      = 1e-6#np.float(raw_input("dt = "))
tend    = 1.7#np.float(raw_input("final time = "))
f_count = 300#np.float(raw_input("number of frames for movie = "))
#
# initialize
#
s, A = 0.04, 0.5
u0   = A*np.exp(-((x/2./s)**2+(y/2./s)**2))       # gaussian
u0   = A*np.exp(-(((x+0.7)/2./s)**2)) + 0*y       # gaussian
v0   = 0.0*(x+y)
u0[0,:],u0[N+1,:],u0[:,0],u0[:,N+1] = 0.,0.,0.,0. # set B.C.
v0[0,:],v0[N+1,:],v0[:,0],v0[:,N+1] = 0.,0.,0.,0.

#mu   = 0.0*(x+y)+1.== dxmu = 0.0*(x+y)
#dxmu = 2./4./np.cos(x/4.)**2+0.0*y
a0   = -2.
ax,ay= -2*0.05, -2*0.2
mu   = a0*np.exp(x**2/ax+(y-0.7)**2/ay) #+0.5    # gaussian
dxmu = a0*2*x/ax*mu
dymu = a0*2*y/ay*mu
mu   = mu + 2.5  # ensure > 0
#rho  = 0.0*(x+y)+1. 
#
# plot initials if True
#
def int_plot(a):
    if a == True:
        plt.figure()
        plt.subplot(221)
        plt.contourf(x, y, u0,cmap=plt.cm.bone)
        plt.title('int-u');plt.colorbar()
        plt.subplot(222)
        plt.contourf(x, y, mu,cmap=plt.cm.BuGn)
        plt.title(r'int-$\mu$');plt.colorbar()
        plt.subplot(223)
        plt.contourf(x, y, dxmu,cmap=plt.cm.inferno)
        plt.title(r'$\partial_x\mu$');plt.colorbar()
        plt.subplot(224)
        plt.contourf(x, y, dymu,cmap=plt.cm.inferno)
        plt.title(r'$\partial_y\mu$');plt.colorbar()
        plt.tight_layout(); plt.show()
    return
int_plot(True)	
#
def init (u0, v0):
# take an euler-rich step
# returns the field at 2 time steps
    u  = u0.copy()            # intitialize/preserve B.C.
    v  = v0.copy()
	
	# define shifted u's
    uxp   = u0[2:N+2 , 1:N+1] # x_i+1
    uxm   = u0[0:N   , 1:N+1] # x_i-1
    uyp   = u0[1:N+1 , 2:N+2] # y_i+1
    uym   = u0[1:N+1 , 0:N]   # y_i-1
    uu    = u0[1:N+1 , 1:N+1] # x_i, y_i
    #rhoi  = rho[1:N+1, 1:N+1] # p_i, P_i
    mui   = mu[1:N+1, 1:N+1]  # mu_i, mu_i
    dxmui = dxmu[1:N+1,1:N+1] # dxmu_i, dxmu_i
    dymui = dymu[1:N+1,1:N+1] # dymu_i, dymu_i

    # coefficiants
    c1 = dxmui/(2*ds)#rhoi*ds)
    c2 = dymui/(2*ds)#rhoi*ds)
    c3 = mui/(ds**2)#rhoi*ds**2)
	
    term1 = c1*(uxp - uxm)
    term2 = c2*(uyp - uym)
    term3 = c3*(uxp + uxm + uyp + uym - 4*uu)
    udd   = term1+term2+term3
    v[1:N+1,1:N+1] = v0[1:N+1,1:N+1] + 0.5*udd*dt 
    u[1:N+1,1:N+1] = u0[1:N+1,1:N+1] + v[1:N+1,1:N+1]*dt
    return u0, u
#
def leap(u0,u):
# function to take a leapfrog step
# returns field at two prior steps
    un = u.copy()           # intitialize and preserve B.C.
	
	# define shifted u's
    uxp = u[2:N+2 , 1:N+1]    # x_i+1
    uxm = u[0:N   , 1:N+1]    # x_i-1
    uyp = u[1:N+1 , 2:N+2]    # y_i+1
    uym = u[1:N+1 , 0:N]      # y_i-1
    uu  = u[1:N+1 , 1:N+1]    # x_i, y_i
    uu0 = u0[1:N+1 , 1:N+1]   # x0_i, y0_i
    #rhoi= rho[1:N+1, 1:N+1]   # p_i, p_i
    mui   = mu[1:N+1, 1:N+1]  # mu_i, mu_i
    dxmui = dxmu[1:N+1,1:N+1] # dxmu_i, dxmu_i
    dymui = dymu[1:N+1,1:N+1] # dymu_i, dymu_i
	
	# coefficiants
    c1 = dt**2*dxmui/(2.*ds)#rhoi*ds)
    c2 = dt**2*dymui/(2.*ds)#rhoi*ds)
    c3 = dt**2*mui/(ds**2)#rhoi*ds**2)
	
    term1 = c1*(uxp - uxm)
    term2 = c2*(uyp - uym)
    term3 = c3*(uxp + uxm + uyp + uym -4*uu)
    term4 = 2*uu - uu0
    un[1:N+1,1:N+1] = term1+term2+term3+term4
    return u, un
#
# loop thru time
#
files  = []; i = 0
levels = np.linspace(-0.9*A,0.9*A,20)          # contour levels
T   = int(tend/dt)                     # number of iterations 
mod = int(T/float(f_count))            # modulus --> saving frame
print 'mod = ',mod,', T = ',T
#
u0, u = init(u0,v0)
for t in range(T):
	# plotting for movie
    if t % mod == 0:
        plt.contourf(x,y,u0,levels,cmap=cmap,extend='both')
        plt.colorbar()
        fname = '_tmp%03d.png' % i; i += 1
        print 'Saving frame', fname
        plt.savefig(fname);files.append(fname);plt.clf()
    u0, u = leap(u0, u)
#
# make movie
#
import os
print('Making movie animation.mpg - this may take a while ')
os.system("ffmpeg -i _tmp%03d.png -c:v libx264 -r 10 -pix_fmt yuv420p dt6_17.mp4")

# cleanup
for fname in files:
    os.remove(fname)
