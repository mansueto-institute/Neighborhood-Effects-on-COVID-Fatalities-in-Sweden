#import pylab
import csv
from matplotlib import *
from pylab import *
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from colorsys import hsv_to_rgb
from scipy.interpolate import make_interp_spline


def linreg(X, Y):
    """
        Summary
        Linear regression of y = ax + b
        Usage
        real, real, real = linreg(list, list)
        Returns coefficients to the regression line "y=ax+b" from x[] and y[], and R^2 Value
        """
    if len(X) != len(Y):  raise ValueError("unequal length")
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in zip(X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
    meanerror = residual = 0.0
    for x, y in zip(X, Y):
        meanerror = meanerror + (y - Sy/N)**2
        residual = residual + (y - a * x - b)**2
    RR = 1 - residual/meanerror
    ss = residual / (N-2)
    Var_a, Var_b = ss * N / det, ss * Sxx / det
    return a, b, RR, Var_a, Var_b


filename='../data/Figure_S9_to_S12.csv' 
f=open(filename, 'r')
reader=csv.reader(f,delimiter=',')

dead_covid=[]
dead=[]
pop=[]
name=[]


for row in reader:
	if (row[2].isdigit()):
		if (float(row[2])==0.):
			print(row[0],row[1],row[2],row[3],row[9])
		name.append(row[0])
		dead_covid.append(float(row[2]) )
		dead.append( float(row[1]) )
		pop.append(float(row[9]))

f.close()


################### first figure: Scaling of Deaths with population (neighborhoods)

print ('')
logpop=[]
logdead=[]
nn=[]
for i in range(len(pop)):
	#if (dead[i]==0.):
	#	print(name[i])
	if (dead[i] >0. ):
		logpop.append( np.log10(pop[i]) )
		logdead.append( np.log10(dead[i]) )
		nn.append(name[i])

xx=logpop
yy=logdead

plt.plot(xx,yy,'o',c='#0066ff',markeredgecolor='white',ms=5,lw=2,label='average margin',alpha=0.4)

gradient, intercept, r_value, var_gr, var_it = linreg(xx,yy)
print('Scaling of deaths, OLS fit parameters, Total:')
print('exponent= ',gradient,", 95 % CI = [",gradient- 2.*np.sqrt(var_gr),",",gradient+2.*np.sqrt(var_gr),"]")
print ('log10 intercept= ',intercept,", 95 % CI = [",(intercept- 2.*np.sqrt(var_it)),",",(intercept+2.*np.sqrt(var_it)),"]")


fitx=np.arange( min(xx)-0.1,max(xx)+0.1,0.1,dtype=float )
fity=intercept + fitx*gradient
plt.plot(fitx,fity,'-', c='blue', linewidth=2, alpha=0.7,label='OLS scaling best fit')

xlabel('neighborhood population',fontsize=20)
ylabel('total deaths',fontsize=20)
#plt.legend(loc='upper left')
plt.tight_layout()
str='Figure_neigh_scaling_deaths.pdf'
savefig(str,dpi=300)
plt.show()

###################### Figure_S9  ###################################

plt.clf()
print ('')
#for i in range(len(pop)):
#	if (dead[i]==0.):
#		print(name[i])

yy=dead_covid
xx=dead

plt.plot(xx,yy,'o',c='#0066ff',markeredgecolor='white',ms=6,lw=2,label='average margin',alpha=0.5)

gradient, intercept, r_value, var_gr, var_it = linreg(xx,yy)
print('COVID deaths vs deaths; OLS fit parameters, Total:')
print('gradient= ',gradient,", 95 % CI = [",gradient- 2.*np.sqrt(var_gr),",",gradient+2.*np.sqrt(var_gr),"]")
print ('intercept= ',intercept,", 95 % CI = [",(intercept- 2.*np.sqrt(var_it)),",",(intercept+2.*np.sqrt(var_it)),"]")


fitx=np.arange( min(xx)-0.1,max(xx)+0.1,0.1,dtype=float )
fity=intercept + fitx*gradient
plt.plot(fitx,fity,'-', c='blue', linewidth=2, alpha=0.9,label='OLS scaling best fit')


#plt.ylim(0.01,0.015)
ylabel('COVID deaths',fontsize=20)
xlabel('total deaths',fontsize=20)
#plt.legend(loc='upper left')
plt.tight_layout()
str='Figure_S9.pdf'
savefig(str,dpi=300)
plt.show()


#########  Figure_S10 ###################


plt.clf()
print ('')

covidpp=[]
deadpp=[]
for i in range(len(pop)):
		covidpp.append( dead_covid[i]/pop[i] )
		deadpp.append( dead[i]/pop[i] )

xx=deadpp
yy=covidpp

plt.plot(xx,yy,'o',c='#0066ff',markeredgecolor='white',ms=6,lw=2,label='average margin',alpha=0.6)

gradient, intercept, r_value, var_gr, var_it = linreg(xx,yy)
print('COVID deaths vs deaths; OLS fit parameters, Total:')
print('gradient= ',gradient,", 95 % CI = [",gradient- 2.*np.sqrt(var_gr),",",gradient+2.*np.sqrt(var_gr),"]")
print ('intercept= ',intercept,", 95 % CI = [",(intercept- 2.*np.sqrt(var_it)),",",(intercept+2.*np.sqrt(var_it)),"]")

fitx=np.arange( min(xx)-0.001,max(xx)+0.001,0.1,dtype=float )
fity=intercept + fitx*gradient
plt.plot(fitx,fity,'-', c='blue', linewidth=2, alpha=0.7,label='OLS scaling best fit')


#plt.ylim(0.01,0.015)
ylabel('COVID deaths/person',fontsize=20)
xlabel('total deaths/person',fontsize=20)
#plt.legend(loc='upper left')
plt.tight_layout()
str='Figure_S10.pdf'
savefig(str,dpi=300)
plt.show()



###################### Figure_S11  ###################################

plt.clf()
print ('')

logcovid=[]
logdead=[]
for i in range(len(pop)):
	#if (dead[i]==0.):
	#	print (name[i],dead[i],dead_covid[i],pop[i])
	if (dead_covid[i] >0.):
		logcovid.append( np.log10( dead_covid[i] ) )
		logdead.append( np.log10( dead[i] ) )

xx=logdead
yy=logcovid

plt.plot(xx,yy,'o',c='#0066ff',markeredgecolor='white',ms=6,lw=2,label='average margin',alpha=0.3)

gradient, intercept, r_value, var_gr, var_it = linreg(xx,yy)
print('log covid deaths vs log deaths; OLS fit parameters, Total:')
print('exponent= ',gradient,", 95 % CI = [",gradient- 2.*np.sqrt(var_gr),",",gradient+2.*np.sqrt(var_gr),"]")
print ('log intercept= ',intercept,", 95 % CI = [",(intercept- 2.*np.sqrt(var_it)),",",(intercept+2.*np.sqrt(var_it)),"]")

#fitx=np.arange( min(xx)-0.1,max(xx)+0.1,0.1,dtype=float )
#fity=intercept + fitx*gradient
#plt.plot(fitx,fity,'-', c='blue', linewidth=2, alpha=0.7,label='OLS scaling best fit')


#plt.ylim(0.01,0.015)
ylabel('COVID deaths/person',fontsize=20)
xlabel('total deaths/person',fontsize=20)
#plt.legend(loc='upper left')
plt.tight_layout()
str='Figure_S11.pdf'
savefig(str,dpi=300)
plt.show()


###################### Figure_S12  ###################################

plt.clf()
print ('')

logcovidpp=[]
logdeadpp=[]
for i in range(len(pop)):
	#if (dead[i]==0.):
	#	print (name[i],dead[i],dead_covid[i],pop[i])
	if (dead_covid[i] >0.):
		logcovidpp.append( np.log10( dead_covid[i]/pop[i] ) )
		logdeadpp.append( np.log10( dead[i]/pop[i] ) )

xx=logdeadpp
yy=logcovidpp

plt.plot(xx,yy,'o',c='#0066ff',markeredgecolor='white',ms=6,lw=2,label='average margin',alpha=0.6)

gradient, intercept, r_value, var_gr, var_it = linreg(xx,yy)
print('COVID deaths vs deaths; OLS fit parameters, Total:')
print('gradient= ',gradient,", 95 % CI = [",gradient- 2.*np.sqrt(var_gr),",",gradient+2.*np.sqrt(var_gr),"]")
print ('intercept= ',intercept,", 95 % CI = [",(intercept- 2.*np.sqrt(var_it)),",",(intercept+2.*np.sqrt(var_it)),"]")

#fitx=np.arange( min(xx)-0.001,max(xx)+0.001,0.1,dtype=float )
#fity=intercept + fitx*gradient
#plt.plot(fitx,fity,'-', c='blue', linewidth=2, alpha=0.7,label='OLS scaling best fit')


#plt.ylim(0.01,0.015)
ylabel('COVID deaths/person',fontsize=20)
xlabel('total deaths/person',fontsize=20)
#plt.legend(loc='upper left')
plt.tight_layout()
str='Figure_S12.pdf'
savefig(str,dpi=300)
plt.show()






