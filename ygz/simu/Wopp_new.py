import aipy as a, numpy as np, capo as C, pylab as plt
from scipy import signal
import seaborn as sns
import w_opp
from joblib import Parallel, delayed
sns.set_context("paper", font_scale=2)
sns.set(style="ticks", color_codes=False,font='DejaVu Serif', font_scale=1.5)
plt.rc('axes', linewidth=1.5)


N_universe = 1
bl1 = (96, 102)
bl2 = (96, 94)
T0 = 2455700.5
fqs = np.linspace(0.145, 0.155, 10)
T1 = np.arange(2455700.49,2455700.51,0.0002)
global WS
WS = w_opp.OppSolver(fqs=fqs, T1=T1)
settings = [(bl1,bl1,True, 0), (bl1, bl1,False, 0),
			(bl1,bl2,True, 'auto'), (bl1,bl2,False, 'auto')]

def run_opp(s):
	print s
	return WS.opp(bl1=s[0], 
				bl2=s[1], 
				sky=s[2], 
				rephase=s[3], 
				delay=True, 
				return_series=True)

res = Parallel(n_jobs=len(settings))(delayed(run_opp)(s) for s in settings)

res_vis0 = res[0]
res_the0 = res[1]
norm_vis = np.amax(np.abs(res_vis0))
norm_the = np.amax(np.abs(res_the0))


res_vis = res[2]
res_the = res[3]
res_vis0 /= norm_vis
res_the0 /= norm_the
res_vis /= norm_vis
res_the /= norm_the
plt.rc('text', usetex=False)
f,axes = plt.subplots(2,1, sharex=True)
axes[0].plot(T1-T0, res_vis0.real, label=r'$<V*V>$')
axes[0].plot(T1-T0, res_the0.real, label=r'$\Theta$')
axes[0].legend()
axes[1].plot(T1-T0, res_vis.real, label=r'$<V*V>$')
axes[1].plot(T1-T0, res_the.real, label=r'$\Theta$')
axes[1].legend()
plt.setp(axes[0].get_xticklabels(), visible=False)
axes[1].set_xlabel('Offset (Sidereal Days)')
plt.show()
import IPython; IPython.embed()