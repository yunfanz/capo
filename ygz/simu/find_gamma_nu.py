from w_opp import *
import itertools
from scipy.signal import argrelextrema
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
from skimage.restoration import unwrap_phase
import seaborn as sns
from matplotlib.font_manager import FontProperties
import aipy as ap
sns.set_context("paper", font_scale=2)
sns.set(style="ticks", color_codes=False,font='DejaVu Serif', font_scale=1.5)
plt.rc('axes', linewidth=1.5)
font = FontProperties()

fqs = np.arange(0.145, 0.156, 0.001)
T1 = np.arange(2456681.49, 2456681.51, 0.0002)
bl1 = (96,102)
bl2 = (96,94)
# T1 = np.arange(2456681.46, 2456681.54, 0.001)
# bl1 = (0,26)
# bl2 = (0,38)
#WS = OppSolver(T1=T1, fqs=np.array([.15]))
#maxres, T1ac = WS.opp(bl1=(0,26), bl2=(0,38), rephase=0,delay=False, return_series=False)
#T1ac += 2456681.5
#Tbox = np.arange(T1ac-0.1, T1ac+0.1, 0.001)
def get_zenith_dtau(dT):
	cal='psa6622_v003'
	
	aa = ap.cal.get_aa(cal, fqs)
	t0 = 2456681.5
	t1 = t0 + dT
	aa.set_jultime(t0)
	bl1v = aa.get_baseline(bl1[0],bl1[1],'z')
	bl2v = aa.get_baseline(bl2[0],bl2[1],'z')
	top0 = (0,0,1)
	m0 = np.linalg.inv(aa.eq2top_m)
	eq = np.dot(m0, top0)
	bl1proj0 = np.dot(bl1v, top0)
	bl2proj0 = np.dot(bl2v, top0)

	aa.set_jultime(t1)
	top1 = np.dot(aa.eq2top_m, eq)
	bl1proj1 = np.dot(bl1v, top1)
	bl2proj1 = np.dot(bl2v, top1)

	return -(bl1proj1-bl1proj0)


coll = []
def get_res_for_freq(fq, bl1=bl1, bl2=bl2, tracking=False):
	WS = OppSolver(T1=T1, fqs=np.array([fq]), tracking=tracking)
	res = WS.opp(bl1=bl1, bl2=bl2, rephase=0,delay=False, sky=False, return_series=True)
	return res
fig = plt.figure()
#frame1 = fig.add_axes((.1,.3,.8,.6))
for i,tracking in enumerate([True, False]):
	res = Parallel(n_jobs=fqs.size)(delayed(get_res_for_freq)(fq, tracking=tracking) for fq in fqs)
	phi_nu = [np.angle(series[np.argmax(np.abs(series))])[0] for series in res]
	phi_nu = unwrap_phase(np.array(phi_nu))
	coll.append(phi_nu)
	phi_nu -= phi_nu[5]
	if i == 1:
		dT = T1[np.argmax(np.abs(series))] - 2456681.5
		dtau = get_zenith_dtau(dT)
		coll.append(np.pi*2*fqs*dtau)
		print 'Delta Tau=',  (phi_nu[-1]-phi_nu[0])/(fqs[-1]-fqs[0])/2/np.pi, dtau
	else:
		coll.append(np.zeros_like(fqs))
	
	#plt.plot(fqs,phi_nu, label=['tracking', 'drift-scan'][i])

plt.plot(fqs, coll[0]-coll[0][5], c='b', label='Tracking')
plt.plot(fqs, coll[1]-coll[1][5], '--', c='b')
plt.plot(fqs, coll[2]-coll[2][5], c='r', label='Drift-scan')
plt.plot(fqs, coll[3]-coll[3][5], '--', c='r')

#plt.grid()
plt.legend(loc=2)
plt.xlabel('Frequency '+r'$\nu$'+' [GHz]')
plt.ylabel('Phase of '+r'$\Theta_\nu$'+' at\n'+'optimal '+r'$\Delta t$')

# frame2 = fig.add_axes((.1,.1,.8,.2))
# plt.plot(fqs, (coll[2]-coll[2][5])-(coll[3]-coll[3][5]), '--', c='g')
# plt.grid()
plt.show()

# for pair in itertools.combinations(res, 2):
# 	img = pair[0][:,np.newaxis] + pair[1][np.newaxis,:]
# 	print np.unravel_index(np.argmax(np.abs(img)), img.shape)

import IPython; IPython.embed()
