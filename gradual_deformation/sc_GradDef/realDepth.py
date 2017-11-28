#!/usr/bin/python
#
"""
This script converts a set of 2D velocity model realizations generated 
via sequential gaussian simulation with Stanford Geostatistical 
Modeling Software (SGeMS) to a .npy. In short these realizations share 
the same spatial variogram (slowness). 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as op

#%% Functions

def meshLookup(Mesh, x_trial, z_trial):
    """This function takes the mesh properties and a (x,y) location and returns 
    the index of the cell that (x,y) belongs to
    
    input: Mesh is a 4x1 list [x cell size, z cell size, num of cells in x, num 
        of cells in z]
    input: x_trial is the x position (in [m]) that want to find the cell index 
    input: z_trial is the z position (in [m]) that want to find the cell index 
    output: row is the matrix index of z_trial
    output: col is the matrix index of x_trial
    
    example:
    row, col = meshLookup(np.array([25,25,60,80]), 250, 250)
    and
    Matrix[row,col] retrieves the entry at x,y = (250,250) 
    """
    x_spacing = Mesh[0]
    z_spacing = Mesh[1]
    xnum = np.int(Mesh[2])
    znum = np.int(Mesh[3])
    col = 999
    row = 999
    for i in range(xnum):
        if x_spacing*i < x_trial <= x_spacing*(i+1):
            col = i
            break
    for j in range(znum):
        if z_spacing*j <= z_trial <= z_spacing*(j+1):
            row = j
            break
    if col == 999:
        print(x_trial)
        print(z_trial)
    if row == 999:
        print(x_trial)
        print(z_trial)
    return row, col
    
def stepAlongPath( xo, zo, xf, zf, stepLength):
    """This function takes a Tx location (xo, zo) and a Rx location (xf, zf)
    and a step length and calculates all positions along the path with 
    equal seperation of stepLength
    
    input: xo is the x location of transmitter
    input: zo is the z location of transmitter
    input: xf is the x location of receiver
    input: zf is the z location of receiver
    input: stepLength is the discretization length along the ray path
    
    output: nupmy array, L x 2, of x,y points equidistant along ray path
    
    example: 
    path = stepAlongPath(0,1000,1500,0,2)
    """
    assert xf - xo > 0
    assert stepLength > 0
    
    num_steps = np.linalg.norm([xf-xo,zf-zo])/stepLength
    delx = (xf-xo)/num_steps
    delz = (zf-zo)/num_steps
    
    xtrial = xo
    ztrial = zo
    path = []
    if delz > 0:
        while np.linalg.norm([xf-xtrial,zf-ztrial]) > stepLength and xtrial < xf and ztrial < zf:
            xtrial += delx
            ztrial += delz
            path.append([xtrial,ztrial])
    elif delz < 0:
        while np.linalg.norm([xf-xtrial,zf-ztrial]) > stepLength and xtrial < xf and ztrial > zf and ztrial + delz >0:
            xtrial += delx
            ztrial += delz
            path.append([xtrial,ztrial])  
            if ztrial< 0:
                print(xo)
                print(zo)
                print(xf)
                print(zf)
    else:   # delz  =0
        while np.linalg.norm([xf-xtrial,zf-ztrial]) > stepLength and xtrial < xf:
            xtrial += delx
            ztrial += delz
            path.append([xtrial,ztrial])  
    
    return np.array(path)
    
def stepAlongMesh(path, Mesh,stepLength):
    """ This function takes an n x 1 array of path points and returns the mesh indices
    associated with those points
    
    input: path is a L x 2 numpy array, is the output of stepAlongPath
    input: Mesh is a 4x1 list [x cell size, z cell size, num of cells in x, num 
        of cells in z]
    input: stepLength is the discretization length along the ray path
    
    output: Sens is the row entry of a Sensitivity matrix for this ray
    
    example:
    SensRow = stepAlongMesh(path,mesh,stepLength)
    """
    xnum = np.int(Mesh[2])
    znum = np.int(Mesh[3])    
    
    n = np.shape(path)[0]
    Sens = np.zeros((znum,xnum))
    for i in range(n):
        row, col = meshLookup(Mesh, path[i,0], path[i,1] )
        Sens[row,col]+= stepLength
    
    return Sens
  
def mat2vec(mat):
    """ reshape a matrix to a vector by convention of bottom up, inverse 
    operation of vec2mat
    input: matrix [number z cells x number of x cells[]
    output: M*N x 1 vector 
    
    example:
    v = mat2vec(B)
    """
    vertical_dim = np.shape(mat)[0]
    horizontal_dim = np.shape(mat)[1]
    vec = np.zeros(vertical_dim*horizontal_dim)
    for i in range(vertical_dim):
        # iterate over z
        for j in range(horizontal_dim):
        # iterate over x
            vec[i*horizontal_dim+j] = mat[vertical_dim-i-1,j]
    return vec
    
def vec2mat(vec):
    """ This function reshapes vector to array by bottom up convention, inverse
    operation of mat2vec
    input: vec is [N*M x 1 array] 
    input: xnum is the number of cells in x direction
    input: znum is the number of cells in z direction
    
    output: mat is the matrix [znum x xnum] representation of model vector
    
    example:
    B2 = vec2mat(v, 60, 80)
    """
    mat = np.zeros((znum,xnum))
    for i in range(znum):
        #iterate over depths
        for j in range(xnum):
            #iterate over x's
            mat[znum-i-1,j]=vec[i*xnum+j]
    return mat
    
def forwModel(survey, mesh, stepLength):
    """This functions peforms the forward modelling for a mesh, model and survey
    and return the travel time for each ray
    
    input: survey is [N x 4] array of transmitter locations and receiver location,
    each row is a single ray
    input: mesh is a 4x1 list [x cell size, z cell size, num of cells in x, num 
        of cells in z]
    input: stepLength is the discretization length along the ray path
    input: model vector of SLOWNESS, [M x 1]
    
    output: sens_row is the linear forward operator, B
    
    example:
    sens1 = forwModel(Survey1, mesh, 1)
    
    """
    
    num_rays = np.shape(survey)[0]
    sens_row = np.zeros((num_rays, np.int(mesh[2]*mesh[3])))
    traveltimes = np.zeros(num_rays)
    for i in range(num_rays):
        path_temp = stepAlongPath(survey[i,0],survey[i,1],survey[i,2],survey[i,3],stepLength)
        sens_temp = stepAlongMesh(path_temp, mesh, stepLength)
        sens_row[i,:] = mat2vec(sens_temp)
    return sens_row
  
def readReal(datFpath):
	""" Function that reads in raw SGeMS output
	creating an array xnum*znum X realizations"""
	with open(datFpath,"r") as ins:
		lines=[]
		for line in ins:
			lines.append(line)
	numR=len(np.array(lines[-1].split(),dtype='float64'))
	numL=len(lines)-(numR+2)
	R_int=np.zeros([numL,numR])
	for i in range((numR+2),len(lines)):
		tmp=lines[i].split()
		tmp=np.array(tmp,dtype='float64')
		jj=i-(numR+2)
		R_int[jj,:]=tmp
	return R_int

def transformReal(r1):
	""" Function that takes realization array created with readReal()
    	and converts to x-space"""

	r1=np.fliplr(vec2mat(r1))
	r1=mat2vec(r1)
	r1=r1*sigmaWell+muWell
	return r1	
	
def mkplot(gaussHistory,i):
	""" Function used to plot intermediate slowness fields"""

	Sfield=vec2mat(gaussHistory)
	Sstacked=np.insert(Sfield,-1,values=wellData,axis=1)
	
	plt.figure(figsize=(6,6))
	plt.imshow(Sstacked,cmap=plt.cm.seismic,interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax);
	plt.colorbar();plt.xlabel('x cell');plt.ylabel('z cell')
	plt.title('Model- it.'+str(i))
	fname=pre+'\images\model'+str(i+1)
	plt.savefig(fname)
	plt.close();print 'fig-',str(i+1),'-saved'

def objectiveFunc(rD):
	# create the trial
	r_temp=r0*np.cos(np.pi*rD)+r_trial*np.sin(np.pi*rD)
	# forw model the trial
	m_trial=r_temp
	
	t_trial=Sensitivity.dot(m_trial)
	# calc the misfit
	mis_trial=np.linalg.norm(ttime-t_trial)
	
	return mis_trial	
		
#%% import and define data

# grid params from MATLAB 
z=range(150,351);t=range(0,1501) # depth (m), and time (ms)
slowness=0.00033 # reference speed from MATLAB picker
xwell=1.206757668289048e+02
ywell=25.405600113841800;zwell=429

# load well data... these are slowness values with respect to depth interpolated from lithology logs
# slowness values were used as conditioning data in the statistical slowness models
wellData=np.genfromtxt('..\sc_SGeMS\DatFileFinished\mkSGeMS\Well56Aslow-ABV.dat',unpack=True,skip_header=6,dtype='float64')[-1]
muWell=np.mean(wellData);sigmaWell=np.std(wellData)
vmin,vmax=muWell-2*sigmaWell,muWell+2*sigmaWell

# import all shot positions
shot,xs,ys,zs=np.genfromtxt('shotPosition.txt',unpack=True,skip_header=1)

iii=3 # choose shot (3=shot30124)
xshot=xs[iii];yshot=ys[iii];zshot=zs[iii]
shotStr=str(shot[iii])[:5]

# import firt arrival times picked in MATLAB
tshot=np.genfromtxt('ttime'+shotStr+'.txt',skip_header=1)

# set up grid
xydis=np.array([xwell-xshot,ywell-yshot]);xydis=np.sqrt(sum(xydis**2))
print '  shot= ',shot[iii], ', mapDis=',xydis
xsize,zsize,xnum,znum=xydis/199.,1.021,199,350
mesh=np.array([xsize,zsize,xnum,znum])	

# load sensitivity and realizations
""" The sensitivitie.npy were created in 
another script. It is a xnum*znum length vector with 
values:=ray distance traveled thru square (i.e. multiply
the sensitivity matrix with slowness field to forward model
for travel times).
The realizations are created in SGeMS and converted to x-space 
in the script svReal, and saved to realization.npy """

sens=np.load('Sensitivities\sens124.npy') # sensitivity matrix
Rpath='..\sc_SGeMS\Realizations\saved\R_124_00300.npy'
R=np.load(Rpath) # set of realizations

#%% peform gradual deformation!

Sensitivity=np.copy(sens)
ttime=np.copy(tshot)

# forward initial model 
r_original=R[:,0] # first realization
r0=np.copy(r_original)
t0=Sensitivity.dot(r0) # forward
# Calculate the initial misfit
mis0=np.linalg.norm(ttime-t0)

# intitialize
fval_best=100
numCells=xnum*znum
num_r=np.shape(R)[1]
rD_history=np.zeros(num_r);rD_history[0]=0
misfit_history=np.zeros(num_r);misfit_history[0]=mis0
gauss_history=np.zeros((numCells,num_r));gauss_history[:,0]=r0

print 'GradDef begin shot-',shotStr

for i in range(1,num_r):
    print(str(i) + ': ')
    # generate the next realization to test
    global r_trial
    r_trial = R[:,i]
    
    # call brent optimization
    rD = op.fminbound(objectiveFunc, -1, 1, xtol = 1e-6, maxfun = 100)
    
    # test for improvement
    fval = objectiveFunc(rD)
    if fval < fval_best:
        # new solution is better
        global r0
        r0 = r0*np.cos(np.pi*rD) + r_trial*np.sin(np.pi*rD)
        gauss_history[:,i]=r0
        misfit_history[i]=fval
        fval_best = fval
        print(fval)
    else:
        # solution is not improved
        global r0
        r0 = r0
        misfit_history[i]=fval_best
        gauss_history[:,i]=r0
    # save some more stuff
    rD_history[i]= rD

print 'GradDef complete-',shotStr

#%% plot	

# save plots see all in ./images30124
savePLOT=True
if savePLOT:
	pre='images'+str(shot[iii])[:5]	# for saving images
	for i in range(num_r):                   # plot all optimized models
		S=gauss_history[:,i]
		mkplot(S,i)

	sens_sum = np.zeros(xnum*znum)           # plot sensitivity
	for i in range(len(ttime)):
		sens_sum=sens_sum+Sensitivity[i,:]

	plt.matshow(vec2mat(sens_sum), cmap=plt.cm.hot,vmin=0,vmax=1.4)
	plt.xlabel('x [m]');plt.ylabel('y [m]');plt.title('Rays')
	plt.colorbar()	
	fname=pre+'\sensitivity'
	plt.savefig(fname)
	plt.close()
		
	# misfit
	plt.figure(figsize=(6,6))                 # plot misfit
	plt.plot(misfit_history, label='Misfit')
	plt.xlabel('Iteration');plt.ylabel('L2 Misfit');plt.title('Misfit')
	plt.tight_layout()
	fname=pre+'\misFit'
	plt.savefig(fname)
	plt.close()

	# statistics
	# point by point

	# entire field
	stats=[[],[],[],[]] # mean, std, min, max

	for i in range(num_r):
		S=gauss_history[:,i]
		Sfield=vec2mat(gauss_history[:,i])
		stats[0].append(np.mean(Sfield))
		stats[1].append(np.std(Sfield))
		stats[2].append(np.min(Sfield))
		stats[3].append(np.max(Sfield))
		
	meanS,stdS,minS,maxS=stats[0],stats[1],stats[2],stats[3]

	plt.figure(figsize=(6,6))
	numPlots=4;fontsize=9
	plt.subplot(numPlots,1,1)
	plt.plot(meanS,'gx-',label='arithmetic mean')
	plt.title('statistical slowness-metrics ')
	plt.ylabel('[s/m]')
	plt.legend(fontsize=fontsize)

	plt.subplot(numPlots,1,2)
	plt.plot(stdS,'kx-',label='standard deviation')
	plt.ylabel('[s/m]')
	plt.legend(fontsize=fontsize)


	plt.subplot(numPlots,1,3)
	plt.plot(minS,'cx-',label='minimum')
	plt.ylabel('[s/m]')
	plt.legend(fontsize=fontsize)


	plt.subplot(numPlots,1,4)
	plt.plot(maxS,'rx-',label='maximum')
	plt.ylabel('[s/m]')
	plt.xlabel('iteration')
	plt.tight_layout()
	plt.legend(fontsize=fontsize)

	fname=pre+'/stats'
	plt.savefig(fname)
	plt.close()



	# get stat's based on index in all Sfields

	allFields=np.zeros([znum,xnum,num_r])
	meanS=np.zeros([znum,xnum])
	medianS=np.zeros([znum,xnum])
	varS=np.zeros([znum,xnum])


	allFields=np.zeros([znum,xnum,num_r])
	for i in range(num_r):
		Sfield=vec2mat(gauss_history[:,i])
		allFields[:,:,i]=Sfield

	meanS=np.zeros([znum,xnum])
	medianS=np.zeros([znum,xnum])
	stdS=np.zeros([znum,xnum])
	for i in range(znum):
		for j in range(xnum):
			tmp=allFields[i,j,:]
			meanS[i,j]=np.mean(tmp)
			medianS[i,j]=np.median(tmp)
			stdS[i,j]=np.std(tmp)

	plt.figure(figsize=(6,6))
	plt.imshow(meanS,cmap=plt.cm.seismic,interpolation='nearest',aspect='auto')#,vmin=vmin,vmax=vmax);
	plt.colorbar();plt.xlabel('x cell');plt.ylabel('z cell')
	plt.title('Mean')	
	fname=pre+'\meanS'
	plt.savefig(fname)
	plt.close()

	plt.figure(figsize=(6,6))
	plt.imshow(medianS,cmap=plt.cm.seismic,interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax);
	plt.colorbar();plt.xlabel('x cell');plt.ylabel('z cell');plt.title('MedianS')	
	fname=pre+'\medianS'
	plt.savefig(fname)
	plt.close()	

	plt.figure(figsize=(6,6))
	plt.imshow(stdS,cmap=plt.cm.seismic,interpolation='nearest',aspect='auto');
	plt.colorbar();plt.xlabel('x cell');plt.ylabel('z cell')
	plt.title('standard deviation')	
	fname=pre+'\stdS'
	plt.savefig(fname)
	plt.close()

