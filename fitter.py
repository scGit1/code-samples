#!/User/bin/python
#
#
# Script to perfom leastsquares fit 
#
#
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
#
#
#
def fitter(p0,x,y,func,errfunc,err):
# the fitter function   
    out = leastsq(errfunc,p0,args=(x,y,func,err),full_output=1)
    pfinal = out[0]
    covar  = out[1]
    mydict = out[2]
    mesg   = out[3]
    ier    = out[4]
    resids = mydict['fvec']
    chisq  = np.sum(resids**2)
    degs_frdm     = len(x)-len(pfinal) 
    reduced_chisq = chisq/degs_frdm
    print 'fitter status: ', ier, '-- aka -- ', mesg
    i = 0
    if covar != None:
        for u in pfinal:
            print 'Param', i+1, ': ',u, ' +/- ', np.sqrt(covar[i,i])
            i = i + 1
        print 'reduced chisq',reduced_chisq
    else:
        print 'fitter failed'
    list = [pfinal,covar,mydict,mesg,ier,resids,chisq,degs_frdm,reduced_chisq]
    return list
#
#
def func(x,p):
# function has the form of desired fit
# the fitted parameters are in the list p
# so for an exponetial fit of the form y = ae^(bx)
    z = p[0]*np.exp(p[1]*x)
    return z
#
#	
def errfunc(p,x,y,func,err):
# returns the error between measured y and func(y)
        return (y-func(x,p))/err	
#
#
# import data		
#x0,y0    = np.genfromtxt('C:\Users\spc_c\Desktop\H_Gamma_slit0.txt',unpack=True)
#
# or create
x0    = np.linspace(0,1)
noise = 4*np.random.uniform(-1,1,len(x0))   # add some noise     
y0    = 2*np.exp(3*x0)+ noise
#
# the error in measured y
err = 2*np.ones(len(x0))
#
# first guess at parameters (this can be finicky)
p0  = [1.2,4.3] # (i.e. fit guess is y = 1.2e^(4.3x) )
#
# fit
fit = fitter(p0,x0,y0,func,errfunc,err)
#
# then pull values you want
params     = fit[0]; a=params[0]; b=params[1]
covariance = fit[1]
residuals  = fit[5]
#
# plot fit against data with errorbars 
x = np.linspace(0,1)
yfit = a*np.exp(b*x)
plt.plot(x,yfit, 'c-')
plt.errorbar(x0,y0,yerr=err,fmt='x', ecolor='g', capthick=2)
plt.show()
