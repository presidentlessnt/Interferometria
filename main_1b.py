#! /usr/bin/python3
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from scipy import odr


dD=.5 #\mu m
E1=np.array([[10,20,30,40,50,60,70,80,90,100,110,120,130],[6,8,11,15,18,20,24,27,31,34,37,40,43]]).T #D,N

D1=np.zeros((len(E1),4))

print(E1[0])

D1[:,0]=E1[:,0]
D1[:,1]=E1[:,1]*2
D1[:,2]=E1[:,0]*0
D1[:,3]=E1[:,0]*0+.5*2

# Lineal function
def func(p,x):
    b,c = p
    return b*x+c

# Model object
quad_model = odr.Model(func)


# Create a RealData object
data = odr.RealData(D1[:,0], D1[:,1], sy=D1[:,3])
#data = odr.RealData(D1[:,0], D1[:,1], sx=D1[:,2], sy=D1[:,3])


# Set up ODR with the model and data.
odr = odr.ODR(data, quad_model, beta0=[630., 1.])

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


slope=popt[0]*1e3
pmslope=perr[0]*1e3

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
errorbar(D1[:,0], D1[:,1], yerr=D1[:,3], hold=True, ecolor='k', fmt='none', label='Data')
xlabel('N')
ylabel('$2 \cdot d$')

plot(x_fit, fit, 'r', lw=2, label='Ajuste lineal \n $\lambda_0:%1.1f \pm %1.1f\;nm$'%(slope,pmslope))
ax.fill_between(x_fit, fit_up, fit_dw, alpha=.25, label='Intervalo 5-sigma')
legend(loc=4)
grid(ls='--',color='grey',lw=.5)
ylim(0,100)
xlim(0,140)
show()

