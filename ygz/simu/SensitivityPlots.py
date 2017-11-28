import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from numpy.random import random

sns.set_context("paper")
#sns.set(style="ticks", color_codes=True,font='DejaVu Serif', font_scale=2)
#plt.rc('axes', linewidth=2.5)

#FILE = 'corr_res.csv'
#FILES = ['HERA_350_pm.csv', 'HERA_243_pm.csv', 'HERA_128_pm.csv', 'HERA_37_pm.csv','PAPER_128_pm.csv']
#FILES = ['HERA_350_all.csv', 'HERA_243_all.csv', 'HERA_128_all.csv', 'HERA_37_all.csv','PAPER_128_all.csv']
#FILES = ['PAPER_128_pm.csv', 'HERA_350_pm.csv', 'HERA_128_pm.csv']
#LABELS = ['HERA350', 'HERA243', 'HERA128', 'HERA37', 'PAPER128']
#LABELS = ['PAPER128','HERA350', 'HERA128']

#FILES = FILES[::-1]; LABELS = LABELS [::-1]
def gen_color(l=1):
	colors = []
	for i in range(l): colors.append((random(),random(),random()))
	return np.array(colors)


def pairplot(Theta_min=0):
	FILES = ['PAPER_128_pm.csv', 'HERA_350_pm.csv', 'HERA_128_pm.csv']
	LABELS = ['PAPER-128','HERA-350', 'HERA-128']

	sns.set(style='ticks', font_scale=1.5,font='DejaVu Serif')
	dflist = []
	for i, file in enumerate(FILES):
		df = pd.read_csv(file)
		df['Array'] = LABELS[i]
		df['$\Theta_{bb}$'] = df['peak']
		df['$\Theta_{bb}$'] /= np.amax(df['peak'])
		df[r'$\bar{b}$'] = df['bl1']
		df['rho0'] = 0.001*40/df['bl1']
		df[r"$\widetilde{\Theta}_{bb^\prime}$"] = np.sqrt(df['mult'])*df['peak']/np.sqrt(1+df['rho0']*2*np.sqrt(df['mult']))
		df[r"$\widetilde{\Theta}_{bb^\prime}$"] /= np.amax(df[r"$\widetilde{\Theta}_{bb^\prime}$"])
		df = df.loc[df[r"$\widetilde{\Theta}_{bb^\prime}$"]>Theta_min]
		df['$\Delta t_{bb^\prime}$'] = df['dT']*24
		dflist.append(df)

	df = pd.concat(dflist)
	plt.locator_params(axis='x', nticks=10)
	g = sns.pairplot(df,hue='Array',vars=['$\Delta t_{bb^\prime}$',r"$\widetilde{\Theta}_{bb^\prime}$",r'$\bar{b}$'],
		plot_kws={"s":30}, diag_kws={'histtype':"step", "linewidth":3})
	for i, j in zip(*np.triu_indices_from(g.axes, 1)):
		g.axes[i, j].set_visible(False)
	#g.axes[0,0].set_visible(False)
	plt.legend(loc=2)

	# g = sns.PairGrid(df,hue='Array',vars=['$dT$','$\Theta$','$\widetilde{\Theta}$','$L$'])
	# g = g.map_diag(plt.hist, histtype="step", linewidth=3)
	# #g = g.map_upper(sns.kdeplot)
	# g = g.map_upper(plt.scatter, alpha=0.5, s=10)
	# g = g.map_lower(plt.scatter, alpha=0.5, s=30)
	# g = g.add_legend()

	ax0 = g.axes[0,0]
	start, end = ax0.get_xlim()
	#start+=0.1999999; end-=0.199999
	ax0.set_xticks(np.linspace(start, end, 3))
	ax0.set_yticklabels([])
	ax1 = g.axes[0,1]
	start, end = ax1.get_xlim()
	#start+=0.2; end-=0.2
	ax1.set_xticks(np.linspace(0, 1, 3))
	ax2 = g.axes[0,2]
	start, end = ax2.get_xlim()
	ax2.set_xticks(np.linspace(0, end, 3))

	#import IPython; IPython.embed()
	plt.gcf().subplots_adjust(right=0.9)
	#g.legend.remove()
	# g = sns.PairGrid(df)
	# g = g.map_diag(sns.kdeplot, lw=3)
	# g = g.map_offdiag(sns.kdeplot, lw=1)

def get_imp(df, Theta_min=0.0):
	dft = df.loc[df['Theta']>Theta_min]
	dfeq = dft.loc[dft['sep']==dft['sep2']]
	dfnq = dft.loc[dft['sep']!=dft['sep2']]
	
	totalsumsq = np.sum(dft['Theta']**2)
	eqsumsq = np.sum(dfeq['Theta']**2)
	totalsum = np.sum(dft['Theta'])
	eqsum = np.sum(dfeq['Theta'])
	totalsens = totalsumsq#/totalsum*np.sqrt(len(dft.index))
	eqsens = eqsumsq#/eqsum*np.sqrt(len(dfeq.index))
	improve = (totalsens-eqsens)/eqsens
	return np.sqrt(totalsens), np.sqrt(eqsens)

def sensplot():
	FILES = ['HERA_350_all_save.csv', 'HERA_243_all.csv', 'HERA_128_all.csv', 'HERA_37_all.csv','PAPER_128_all.csv']
	LABELS = ['HERA-350', 'HERA-243', 'HERA-128', 'HERA-37', 'PAPER-128']
	COLORS = gen_color(len(FILES))
	print "========= Statistics of sensitibity contribution =========="
	sns.set(style="ticks", color_codes=True,font='DejaVu Serif', font_scale=2)

	plt.figure()
	plt.rc('axes', linewidth=2)
	for i, file in enumerate(FILES):
		df = pd.read_csv(file)
		df['rho0'] = 0.001*40/df['bl1']
		df['Theta'] = np.sqrt(df['mult'])*df['peak']/np.sqrt(1+df['rho0']*2*np.sqrt(df['mult']))
		df['Theta'] /= np.amax(df['Theta'])
		#df['Theta'] /= 1000

		TL = np.arange(0,1,0.01)
		 
		imp = np.array([get_imp(df, l) for l in TL])
		totalsensL, eqsensL = imp.T
		S = np.amax(totalsensL)
		totalsensL/=S; eqsensL/=S
		plt.plot(TL, totalsensL, label=LABELS[i], color=COLORS[i], linewidth=3)
		plt.plot(TL, eqsensL,'--', color=COLORS[i], linewidth=3)
	plt.legend(loc=3, fontsize=16)
	plt.xlabel("Effective Correlation Cutoff "+r'$\widetilde{\Theta}_{min}$')
	plt.ylabel("Fraction of Array Sensitivity")
	#plt.ylabel(r'$\rho(\widetilde{\Theta}_{min})/\rho(\widetilde{\Theta}_{min}=0.0)$')
	plt.gcf().subplots_adjust(bottom=0.2)
	#plt.rc('axes', linewidth=2)


if __name__=="__main__":
	#sensplot()
	pairplot(0.01)

	plt.show()
	#import IPython; IPython.embed()

