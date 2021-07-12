#! /usr/bin/python3
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from scipy import odr


t=3.00 #cm
dt=.005 #cm
dP=10 #mBar
E1=np.array([[100,200,300,400,500,600,700],[2,4,6,9,10,12,13]]).T #mBar,N

lamb=0.6332399925575554 # \mu m
dlamb=0.007409711713671636 # \mu m

D1=np.zeros((len(E1),4))


D1[:,0]=E1[:,0]
D1[:,1]=E1[:,1]*lamb/2/t*1e-4
D1[:,2]=E1[:,0]*0+dP #mBar
D1[:,3]=E1[:,1]/2*((dlamb/t*1e-4)**2+(-lamb/t**2*1e-2*dt)**2)**.5

# Lineal function
def func(p,x):
    b,c = p
    return b*x+c

# Model object
quad_model = odr.Model(func)


# Create a RealData object
#data = odr.RealData(D1[:,0], D1[:,1], sy=D1[:,3])
data = odr.RealData(D1[:,0], D1[:,1], sx=D1[:,2], sy=D1[:,3])

print(data)

# Set up ODR with the model and data.
odr = odr.ODR(data, quad_model, beta0=[2., .1])

# Run the regression.
out = odr.run()

#print fit parameters and 1-sigma estimates
popt = out.beta
perr = out.sd_beta
print("fit parameter 1-sigma error")
print("———————————–")

for i in range(len(popt)):
    print(str(popt[i])+"+-"+str(perr[i]))

slope=popt[0]*1e7
pmslope=perr[0]*1e7

# prepare confidence level curves

nstd = 5. # to draw 5-sigma intervals
popt_up = popt + nstd * perr
popt_dw = popt - nstd * perr

x_fit = D1[:,0] #np.linspace(min(x), max(x), 100)
fit = func(popt, x_fit)
fit_up = func(popt_up, x_fit)
fit_dw= func(popt_dw, x_fit)

#plot
fig, ax = plt.subplots(1)
errorbar(D1[:,0], D1[:,1], yerr=D1[:,3], xerr=D1[:,2], hold=True, ecolor='k', fmt='none', label='Data')
xlabel('$\Delta P\;[mbar]$')
ylabel('$\dfrac{N*\lambda}{2*d}$')

plot(x_fit, fit, 'r', lw=2, label='Ajuste lineal \n $\dfrac{\Delta n}{\Delta P}:%1.1f \pm %1.1f\;[\cdot 10^{-7}\,mbar^{-1}]$'%(slope,pmslope))
ax.fill_between(x_fit, fit_up, fit_dw, alpha=.25, label='Intervalo 5-sigma')
legend(loc=4)
grid(ls='--',color='grey',lw=.5)
ylim(0,200e-6)
xlim(50,750)
ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
show()
