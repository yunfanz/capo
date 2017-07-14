import aipy as a, numpy as n, pylab as plt, ephem as e
from mpl_toolkits.mplot3d import Axes3D
#Plot tracks of the entire array as the earth rotates
#aa=a.cal.get_aa('psa6622_v001',n.array([.15]))
import seaborn as sns
#plt.style.use('seaborn-deep')
sns.set(style='ticks', font_scale=1.5,font='DejaVu Serif')
aa = a.cal.get_aa('psa6622_v003',n.array([.15]))
nants = len(aa)
nants = 16
rad2deg=180/n.pi
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15)
#src=a.fit.RadioSpecial("Sun")

#aa.set_jultime(2456240.2)
dt = 0.002
TIME = n.arange(2456249.1,2456249.6, dt)
U, V, W = {}, {}, {}

for time in TIME:
    aa.set_jultime(time)
    src.compute(aa)
    # for i in range(nants):
    #     for j in range(i+1,nants):
    for i, j in [(0,26), (0,38), (0,46)]:
        if src.alt>0:
            u,v,w = aa.gen_uvw(i,j,src=src)
            u,v,w = u.flatten().flatten()[0], v.flatten().flatten()[0], w.flatten().flatten()[0]
            U[(i,j)] = U.get((i,j),[])+[u]
            V[(i,j)] = V.get((i,j),[])+[v]
            W[(i,j)] = W.get((i,j),[])+[w]
            #plt.plot(u,v,'.',ms=2,color = (((j*17)%127)/127., (32-i)*17%127/127., (i+j)/357., 1))
            #print u,v
    
        #plt.plot(-u,-v,'ko')
    

#rs = 10**n.arange(1,2.5,rstep)
#rs = 2**(n.arange(3,8,1) +.5)
#for r in rs:
#    th = n.arange(0, 2*n.pi+.02, .01)
#    x,y = r*n.cos(th), r*n.sin(th)
#    plt.plot(x,y,'r-')
# import seaborn as sns
# sns.set_context("paper")
# sns.set(style="ticks", color_codes=True,font='DejaVu Serif', font_scale=2)
# plt.rc('axes', linewidth=2.5)
# plt.figure()
# #import IPython; IPython.embed()
# for k, u in U.iteritems(): 
#     v = V[k]
#     plt.plot(u,v)
# #plt.xlim(-200,200)
# #plt.ylim(-200,200)
# plt.gcf().subplots_adjust(bottom=0.15, left=0.15)

# plt.grid()
# plt.xlabel('u')
# plt.ylabel('v')
C = {(0,26):'r', (0,38):'b', (0,46):'y'}

fig = plt.figure()
ax = fig.add_subplot(211)
for k, u in U.iteritems(): 
     v = V[k]
     #import IPython; IPython.embed()
     ax.plot(u,v, label='parametric curve', c=C[k])

ax.set_xlabel(r'$u$', fontsize=24)
ax.set_ylabel(r'$v$', fontsize=24)
ax.grid()
ax = fig.add_subplot(212)
#import IPython; IPython.embed()
for k, u in U.iteritems(): 
     v = V[k]
     w = W[k]
     #import IPython; IPython.embed()
     ax.plot(v,w, label='parametric curve', c=C[k])

ax.set_xlabel(r'$v$', fontsize=24)
ax.set_ylabel(r'$w$', fontsize=24)
ax.grid()
# ax = fig.add_subplot(212,projection='3d')
# #import IPython; IPython.embed()
# for k, u in U.iteritems(): 
#      v = V[k]
#      w = W[k]
#      #import IPython; IPython.embed()d
#      ax.plot(u,v,w, label='parametric curve', c=C[k])

# ax.set_xlabel('u')
# ax.set_ylabel('v')
# ax.set_zlabel('w')
plt.tight_layout(h_pad=0.1)
plt.show()
