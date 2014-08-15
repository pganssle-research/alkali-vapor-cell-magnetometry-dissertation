'''
Created on Apr 11, 2013

@author: Paul J Ganssle
'''
import matplotlib
matplotlib.use('Cairo')

from matplotlib import rcParams, rc
from matplotlib.pyplot import *;

import scipy.io;
from numpy import *;

from analyzedata.fileio.importdata import importdata as idata;

bluec='#3953a4';
redc = '#b5130c';
greenc = '#007935';

rcParams['text.usetex']=False
rcParams['mathtext.default'] = 'regular';
rcParams['xtick.major.pad'] = '5'
rcParams['font.family'] = 'Myriad Pro';

rc('font',**{'family':'sans-serif','sans-serif':['Myriad Pro']});

# The Virial coefficients of a few alkalis.
Ts = linspace(100, 250, 150)+273.15;

# Pressures in mtorr, first order corrections.
poff = log10(1/760000.);
Toff = log10(9./5.);

#CsV1 = [8.1636, -7329, -1.5];
#RbV1 = [8.2817, -7706, -1.0];
#KV1 = [8.6853, -8513, -1.2];

#CsP1 = pow(10, CsV1[0]-poff + (CsV1[1]*5./9.)/Ts + CsV1[2]*(Toff + log10(Ts)));
#RbP1 = pow(10, RbV1[0]-poff + (RbV1[1]*5./9.)/Ts + RbV1[2]*(Toff + log10(Ts)));
#KP1 = pow(10, KV1[0]-poff + (KV1[1]*5./9.)/Ts + KV1[2]*(Toff + log10(Ts)));

# Next order correction.
CsV2 = [8.232-poff, -4062, -1.3359];
RbV2 = [8.316-poff, -4275, -1.3102];
KV2 = [8.233-poff, -4693, -1.2403];

CsP2 = pow(10, CsV2[0] + CsV2[1]/Ts + CsV2[2]*log10(Ts));
RbP2 = pow(10, RbV2[0] + RbV2[1]/Ts + RbV2[2]*log10(Ts));
KP2 = pow(10, KV2[0] + KV2[1]/Ts + KV2[2]*log10(Ts));

linew = 1.75;

for ii, backend in enumerate(['cairo', 'Qt4Agg']):
	matplotlib.pyplot.switch_backend(backend);

	nf = figure(num=5, figsize=(16,(9/4)*(16/6.5)), dpi=80);

#	[pcs1] = semilogy(Ts-273.15, CsP1, 'b', color=bluec, lw=linew)
	[pcs2] = semilogy(Ts-273.15, CsP2, 'b', color=bluec, lw=linew)

#	[prb1] = semilogy(Ts-273.15, RbP1, 'r', color=redc, lw=linew)
	[prb2] = semilogy(Ts-273.15, RbP2, 'r', color=redc, lw=linew)

#	[pk1] = semilogy(Ts-273.15, KP1, 'g', color=greenc, lw=linew)
	[pk2] = semilogy(Ts-273.15, KP2, 'g', color=greenc, lw=linew)


	xlabel(u'Temperature (\u00B0C)');
	ylabel(r'Vapor Pressure (mTorr)');

	xlim(Ts[0]-273.15, Ts[-1]-273.15)
	#ylim(10**(log10(KP1[0])-0.25), 25*floor(CsP1[-1]/25.0))

#	leg = legend([pcs1, pcs2, prb1, prb2, pk1, pk2],
#				 ['Cesium (Cs)', ' ', 'Rubidium (Rb)', ' ', 'Potassium (K)', ' '], loc=2);
	leg = legend([pcs2, prb2, pk2],
				 ['Cesium (Cs)', 'Rubidium (Rb)', 'Potassium (K)'], loc=2);

	f = gcf();
	f.set_alpha(0);
	leg.legendPatch.set_facecolor('none')
	ax = gca();
	ax.set_alpha(0);

	ax.xaxis.labelpad = 10;
	ax.xaxis.label.set_fontsize(16);
	ax.xaxis.label.set_fontweight('semibold');

	ax.yaxis.label.set_fontsize(16);
	ax.yaxis.label.set_fontweight('semibold');
	ax.yaxis.labelpad = 3;

	tfsize = 14;

	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(tfsize);

	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(tfsize);

	# ax.title.set_fontsize(17);

	ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%s'))
	ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: '{:g}'.format(x)))

	subplots_adjust(left=0.05,right=0.975, top=0.95, bottom=0.125)

	if ii == 1:
		show();
	else:
		savefig('../magnetometer/AlkaliVaporPressures.pdf', facecolor='none', edgecolor='none', transparent='True', format='pdf');

#savefig('../magnetometer/AlkaliVaporPressures.eps', facecolor='none', edgecolor='none', transparent='True', format='eps');
