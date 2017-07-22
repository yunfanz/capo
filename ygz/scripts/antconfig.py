#plots configuration of antenna array
import aipy as a, numpy as n, pylab as p, ephem as e
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import matplotlib.ticker as tic
import seaborn as sns
from matplotlib import collections  as mc


#plt.style.use('seaborn-deep')
sns.set(style='ticks', font_scale=1.5,font='DejaVu Serif')
plt.rc('axes', linewidth=1.5)
aa = a.cal.get_aa('psa6622_v003',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
nants = 128
rad2deg = 180/n.pi
# ltsec = 299792458.  #meters
# X,Y,Z,I = [],[],[],[]
# for i in range(nants):
#     #print i
#     a = aa.ants[i]
#     pos = a.pos*ltsec*1.E-9       #X,Y,Z position in meters
#     X.append(pos[0]); Y.append(pos[1]); Z.append(pos[2]); I.append(i)
antpos = [aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))]
antpos = n.array(antpos) * a.const.len_ns / 100.
X,Y,Z = antpos[:,0], antpos[:,1], antpos[:,2]
X -= n.average(X)
Y -= n.average(Y)
I = n.arange(nants)

fig = p.figure()
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(211)
p.scatter(X,Y, c='black')
ax.set_aspect(1)
setp( ax.get_xticklabels(), visible=False)
#setp( ax.get_yticklabels(), visible=False)
ax.get_yaxis().set_tick_params(direction='in')
ax.get_xaxis().set_tick_params(direction='in')
ax.grid()
ax = fig.add_subplot(212)
g = 112
Xg, Yg, Ig = X[:g], Y[:g], I[:g]
p.scatter(Xg,Yg, c='black')

lines = [[(X[5], Y[5]), (X[32], Y[32])], 
		[(X[40], Y[40]), (X[14], Y[14])],
		[(X[17], Y[17]), (X[32], Y[32])], 
		[(X[17], Y[17]), (X[55], Y[55])]]
c = ['r','r','y','b']

lc = mc.LineCollection(lines, colors=c, linewidths=2)
ax.add_collection(lc)
#for x,y,i in zip(Xg, Yg,Ig):
#    ax.annotate('%s' %i, xy=(x,y), textcoords='data', fontsize=12) # <--
ax.set_xlabel('East Position [m]')
ax.set_ylabel('North Position [m]')
ax.get_yaxis().set_tick_params(direction='in')
ax.get_xaxis().set_tick_params(direction='in')
ax.set_aspect(5.1)
#ax.grid()

#p.tight_layout()
p.show()
