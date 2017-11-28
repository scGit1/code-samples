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

#%% Functions

def mat2vec(mat):
    """ reshape a matrix to a vector by convention of bottom up, inverse 
    operation of vec2mat
    input: matrix [number z cells x number of x cells[]
    output: M*N x 1 vector 
    
    example:
    v = mat2vector(B)
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
    operation of mat2vector
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

# read in set of realizations generated with SGeMS
datFpath='../sc_SGeMS/Realizations/R_124_00300'
R=readReal(datFpath)

# transform realization to data space and convert to .npy
for i in range(R.shape[1]): 
	R[:,i]=transformReal(R[:,i]) 

savePath=datFpath[:-12]+'/saved/'+datFpath[-11:]
np.save(savePath,R) # save transformed realizations


savePLOT=True # plot if desired
if savePLOT:
	R=np.load(savePath+'.npy')
	r=vec2mat(R[:,89]) # the 89th realization (random sample)
	plt.figure(figsize=(6,6))
	plt.imshow(r,cmap=plt.cm.seismic,interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax)
	plt.colorbar();plt.xlabel('x cell');plt.ylabel('z cell')
	plt.title('realization-')
	plt.savefig(savePath)
	plt.show()


