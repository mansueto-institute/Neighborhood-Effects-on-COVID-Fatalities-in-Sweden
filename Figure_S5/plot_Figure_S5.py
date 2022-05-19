#import pylab
import csv
import scipy.stats
from matplotlib import *
from pylab import *
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from colorsys import hsv_to_rgb
from scipy.interpolate import make_interp_spline
from reliability.Fitters import Fit_Everything
from scipy.stats import lognorm

filename='../data/Figure_S5.csv' 
f=open(filename, 'r')
reader=csv.reader(f,delimiter=',')

pc=[]
pop=[]
for row in reader:
	if (row[0].isdigit()):
		#print(row[0],row[1])
		pc.append(row[0])
		pop.append(float(row[1]) )
f.close()

popn=pop/np.sum(pop) # normalize fraction of neighborhoods

#generate sample of size nn
nn=100000
data=[]
for i in range(len(pc)):
	for j in range(round(nn*popn[i])):
		data.append(float(pc[i]))
#print(nn,len(data))
#gf = Fit_Everything(failures=data, show_histogram_plot=True, show_probability_plot=False, show_PP_plot=False, show_best_distribution_probability_plot=False)

# lognormal parameters 
shape,loc,scale = lognorm.fit(data)
print('lognormal paramerters=',shape,loc,scale)  # s is std of log, loc shifts the origin in x-> x-loc, scale is the log-mean = exp(mu)  


x=np.arange(0, 90, 0.1)
pdf = lognorm.pdf(x, shape, loc, scale)
plt.plot(x, pdf, color='black',lw=3,alpha=0.5,label='lognormal fit')


plt.bar(pc,popn,label='fraction neighborhoods') #,'o-',c='green',markeredgecolor='black',ms=8,lw=2,label='average margin',alpha=0.8)

xlabel('Share neighborhood population born in non-Nordic/EU15 countries',fontsize=12)
plt.xticks(np.arange(0, 100, 10))
ylabel('Fraction neighborhoods',fontsize=20)
plt.legend(loc='upper right')
plt.tight_layout()
str='Figure_S5.pdf'
savefig(str,dpi=300)
plt.show()




