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


norm = colors.Normalize(vmin=0, vmax=20)
sm = cm.ScalarMappable(norm, cmap=cm.Paired)
cnt = 1
mk='o'
edge_color='white'


filename='../data/Figure_2.csv' #data from general population survey USA (annual data are available)
f=open(filename, 'r')
reader=csv.reader(f,delimiter=',')

quant=[]
mg=[]
lerr=[]
gerr=[]
for row in reader:
	if ('%' in row[0]):
		print(row[0],row[1])
		quant.append(row[0].replace(' ','').replace('>',''))
		mg.append(float(row[1]))
		lerr.append(float(row[5]))
		gerr.append(float(row[6]))
f.close()

xx1=[0,1,2,3,4,5,6,7.2]
xx=[0,1,2,3,4,5,6,7]

X_Y_Spline1 = make_interp_spline(xx, lerr)
X_Y_Spline2 = make_interp_spline(xx, gerr)

X_ = np.linspace(0., 7.2, 50)
Y1_ = X_Y_Spline1(X_)
Y2_ = X_Y_Spline2(X_)


#plt.fill_between(xx1,lerr,gerr,color='gray', alpha=0.5)

plt.fill_between(X_,Y1_,Y2_,color='gray', alpha=0.3,label='95% confidence interval')
plt.plot(quant,mg,'o-',c='red',markeredgecolor='black',ms=8,lw=2,label='average margin',alpha=0.8)
#plt.plot(pbrack,ieplus,'k-',alpha=0.4)
#plt.plot(pbrack,ieminus,'k-',alpha=0.4)

plt.ylim(0.01,0.017)
xlabel('Share of neighborhood population with low education',fontsize=14)
plt.xticks(rotation=90)
ylabel('Predictive Margin',fontsize=20)
plt.legend(loc='upper left')
plt.tight_layout()
str='Figure_2.pdf'
savefig(str,dpi=300)
plt.show()




