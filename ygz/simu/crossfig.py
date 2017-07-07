#! /usr/bin/env python
import aipy as a, numpy as np, capo as C, pylab as plt
from scipy import signal
from matplotlib.colors import LogNorm
import numpy.ma as ma
#@plt.ion()
#fqs = np.linspace(.1,.2,203)
fq = .15
bl1, bl2 = (0,26),(0,38)
#cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0., '50_57':0.122}
cuedict = {'26_26':0.,'26_38': 0.03300,'26_46': -0.034, '26_50':0.073557,'13_32':0.030557,'13_14':0.066557,'50_59':0.071557}

dT = 0.033
T0 = 2455700.4
aa = a.cal.get_aa('psa6240_v003', np.array([fq]))
aa.set_jultime(T0)

h = a.healpix.HealpixMap(nside=256)
#h = a.healpix.HealpixMap(nside=64)
#h.set_interpol(False)
#ex,ey,ez = h.px2crd(np.arange(h.map.size), ncrd=3)

bl1x, bl1y, bl1z = aa.get_baseline(bl1[0],bl1[1],'z')
bl2x, bl2y, bl2z = aa.get_baseline(bl2[0],bl2[1],'z')

#fng1=1v

#tx,ty,tz = h.px2crd(np.arange(h.map.size), ncrd=3)
img = a.img.Img(200,.5)
tx2,ty2,tz2 = img.get_top((200,200))
SH = tx2.shape
tx,ty,tz = tx2.flatten(),ty2.flatten(),tz2.flatten()
top = np.array([tx,ty,tz], dtype=tx2.dtype)

bm = aa[0].bm_response((tx,ty,tz),pol='I')[0]**2#/np.abs(tz)   #tz is the Jacobian

#bm = ma.masked_where(tz < 0.001, bm)

#bm /= bm.sum()
#h.map = bm

# C.plot.plot_hmap_ortho(h,mx=0,drng=2)
# plt.colorbar()
# plt.show()


# #Create equatorial coordinates of the first frame T0
top = np.array([tx,ty,tz], dtype=tx.dtype)

m = np.linalg.inv(aa.eq2top_m)
ex,ey,ez = np.dot(m, top)
eq = np.array([ex,ey,ez], dtype=ex.dtype)

T1 = T0+dT
aa.set_jultime(T1)
m = aa.eq2top_m
t2x,t2y,t2z = np.dot(m, eq)

bl2_prj = t2x*bl2x + t2y*bl2y + t2z*bl2z
bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
fng1 = np.exp(-2j*np.pi*bl1_prj*fq)
fng2 = np.exp(-2j*np.pi*bl2_prj*fq)
# #fng2=1
bm2 = aa[0].bm_response((t2x,t2y,t2z),pol='I')[0]**2#/np.abs(tz)#*np.abs(tzsave)
bm2 = ma.masked_where(tz < 0.001, bm2)
# #bm = np.ones_like(tx)
# #bm = np.where(tz > 0, bm, 0)
# bm2 = np.where(tz > 0.001, bm2, 0)
# #import IPython; IPythonp.embed()
# #print bm.sum()
# bm2 /= bm2.sum()
bm_fng2 = bm2 * fng2
bm_fng1 = bm * fng1

#h.map = bm_fng2*bm_fng1.conj()


# h.map = bm_fng1
# plt.figure()
# plt.subplot(131)
# C.plot.plot_hmap_ortho(h,mode="real")
# plt.colorbar()
# h.map = bm_fng2
# plt.subplot(132)
# C.plot.plot_hmap_ortho(h,mode="real")
# plt.colorbar()
# h.map = np.abs(bm_fng2*bm_fng1.conj())
# plt.subplot(133)
# C.plot.plot_hmap_ortho(h,mode="real")
# plt.colorbar()
# plt.show()

bf1 = bm_fng1.reshape((400,400))
bf2 = bm_fng2.reshape((400,400))
bfc = (bm_fng2*bm_fng1.conj()).reshape((400,400))
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
fftfq = np.fft.fftshift(np.fft.fftfreq(400, 2./400))
fq0, fq1 = fftfq[0], fftfq[-1]
im1 = ax1.imshow(bf1.real,extent=[-1,1,-1,1])
im2 = ax2.imshow(bf2.real,extent=[-1,1,-1,1])
im3 = ax3.imshow(np.abs(bfc),extent=[-1,1,-1,1])
Z1 = np.abs(np.fft.fftshift(np.fft.fftn(bfc)))
Z2 = np.abs(np.fft.fftshift(np.abs(np.fft.fftn(bf1)*np.fft.fftn(bf2).conj())))
im4 = ax4.imshow(np.abs(np.fft.fftshift(np.fft.fftn(bf1))), 
	extent=[fq0, fq1, fq0, fq1],
	norm=LogNorm(vmin=Z2.max()*1.e-9, vmax=Z2.max()),
	interpolation='nearest')
im5 = ax5.imshow(np.abs(np.fft.fftshift(np.fft.fftn(bf2))),
	extent=[fq0, fq1, fq0, fq1],
	norm=LogNorm(vmin=Z2.max()*1.e-9, vmax=Z2.max()),
	interpolation='nearest')
im6 = ax6.imshow(Z2, 
	norm=LogNorm(vmin=Z2.max()*1.e-9, vmax=Z2.max()), 
	extent=[fq0, fq1, fq0, fq1])
fig.subplots_adjust(right=0.9)
ax1.set_xlabel('l'); ax2.set_xlabel('l'); ax3.set_xlabel('l')
ax1.set_ylabel('m')
ax4.set_xlabel('u'); ax5.set_xlabel('u'); ax6.set_xlabel('u')
ax4.set_ylabel('v')
for ax in [ax2, ax3, ax5, ax6]:
	plt.setp(ax.get_yticklabels(), visible=False)
cbar_ax3 = fig.add_axes([0.925, 0.57, 0.025, 0.3])
cbar_ax6 = fig.add_axes([0.925, 0.13, 0.025, 0.3])
fig.colorbar(im3, cax=cbar_ax3)
fig.colorbar(im6, cax=cbar_ax6)
plt.show()

#import IPython; IPython.embed()

# #IM = a.img.Img(size=200)
# plt.ion()
# img = a.img.Img(200,.5)
# tx2,ty2,tz2 = img.get_top((200,200))
# SH = tx2.shape
# tx2,ty2,tz2 = tx2.flatten(),ty2.flatten(),tz2.flatten()
# top2 = np.array([tx2,ty2,tz2], dtype=tx2.dtype)
# ex2,ey2,ez2 = np.dot(m, top2)
# h.map = bm_fng1
# bmfng_proj = h[ex2,ey2,ez2]
# bmfng_proj.shape = SH
# plt = plt.imshow(bmfng_proj.imag, vmax=2, vmin=-2,origin='lower')
# plt.show()
# #import IPython; IPythonp.embed()