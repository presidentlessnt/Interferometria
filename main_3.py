#! /usr/bin/python3
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from scipy import odr
import modulo as mtd

t=.56 #cm
dt=.005 #cm
dO=.05*np.pi/180 #deg
#E1=np.array([[3.4,4.5,4.8,5.7,6.3,8.0,8.5,9.2,9.9,10.4,10.9,11.4,11.8,12.2,12.3,13.1,13.4,13.9,14.1,14.5,14.9,15.2,15.5,15.9,16.3,16.7,17.0,17.2,17.0,18.1],[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300]]).T #deg,N
E1=np.array([[3.4,4.5,4.8,5.7,6.3,8.0,8.5,9.2,9.9,10.4,10.9,11.4,11.8,12.2,12.3,13.1,13.4,13.9,14.1,14.5,14.9,15.2,15.5,15.9,16.3,16.7,17.0,17.2,18.0,18.1],[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300]]).T #deg,N


lamb=0.6332399925575554 # \mu m
dlamb=0.007409711713671636 # \mu m

D2=np.zeros((len(E1),4))

E1[:,0]*=np.pi/180

for i in range(len(E1)):
    D2[i,0]=2*t*1e4*(1-np.cos(E1[i,0]))-E1[i,1]*lamb #\mu m
    D2[i,1]=(2*t*1e4-E1[i,1]*lamb)*(1-np.cos(E1[i,0])) #\mu m
    D2[i,2]=((2*(1-np.cos(E1[i,0]))*dt*1e4)**2+(-E1[i,1]*(1-np.cos(E1[i,0]))*dlamb)**2+((2*t*1e4-E1[i,1]*lamb)*np.sin(E1[i,0])*dO)**2)**.5
    D2[i,3]=((2*(1-np.cos(E1[i,0]))*dt*1e4)**2+(-E1[i,1]*dlamb)**2+(2*t*1e4*np.sin(E1[i,0])*dO)**2)**.5


# Lineal function
def func(p,x):
    b,c = p
    return b*x+c

# Model object
quad_model = odr.Model(func)


# Create a RealData object
data = odr.RealData(D2[:,0], D2[:,1], sx=D2[:,2], sy=D2[:,3])

# Set up ODR with the model and data.
odr = odr.ODR(data, quad_model, beta0=[1.5, .1])

# Run the regression.
out = odr.run()

#print fit parameters and 1-sigma estimates
popt = out.beta
perr = out.sd_beta
print("fit parameter 1-sigma error")
print("———————————–")
out.pprint()


for i in range(len(popt)):
    print(str(popt[i])+"+-"+str(perr[i]))

slope=popt[0]
pmslope=perr[0]

# prepare confidence level curves
nstd = 5. # to draw 5-sigma intervals
popt_up = popt + nstd * perr
popt_dw = popt - nstd * perr

x_fit = D2[:,0]
fit = func(popt, x_fit)
fit_up = func(popt_up, x_fit)
fit_dw= func(popt_dw, x_fit)


#plot
fig, ax = plt.subplots(1)
errorbar(D2[:,0], D2[:,1], yerr=D2[:,3], xerr=D2[:,2], hold=True, ecolor='k', fmt='none', label='Data')
ylabel(r'$(2t-N\lambda_0)(1-\cos\theta)\;[\mu m]$')
xlabel(r'$2t(1-\cos\theta)-N\lambda_0\;[\mu m]$')

plot(x_fit, fit, 'r', lw=2, label='Ajuste lineal \n $n_g:%1.2f \pm %1.2f$'%(slope,pmslope))
ax.fill_between(x_fit, fit_up, fit_dw, alpha=.25, label='Intervalo 5-sigma')
legend(loc=4)
grid(ls='--',color='grey',lw=.5)
ylim(0,6e2)
xlim(0,4e2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
ticklabel_format(style='sci', axis='x', scilimits=(0,0), useMathText=True)


plt.figure(2)
plot(E1[:,0],D2[:,1]/D2[:,0],'o-',label='Data')
plot(E1[:,0],E1[:,0]*0+slope,'r',lw=2,label='Ajuste lineal')
ylabel('$n_g$')
xlabel(r'$\theta\;[rad]$')
grid(ls='--',color='grey',lw=.5)
ylim(1.4,2.0)
xlim(.05,.35)
legend()

show()

