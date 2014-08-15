"""
Created on Apr 4, 2013

@author: Omega
"""

import scipy.io;
from numpy import *;

from matplotlib import rcParams, rc
from matplotlib.pyplot import *;

import theme.colors as thc
import theme.figure_settings as tfs

from analyzedata.fileio.importdata import importdata as idata;
from os import path

fname = 'cpmg_refocusing'
floc = 'data_sets'
en=1;
fname = '{}_en_{:02.0f}.mat'.format(path.join(floc, fname), en)

out = idata(fname);

t = out['t'];
ns = out['ns'];
o = out['out'];
oe = out['oe'];

tfs.setup_environment()

tau = 2./ns[2:];
oe2e = oe[2:, :];
linewidth=2;
markeredge = 1.5;
markersize = 6;

fsize=(8, 6)

fnum = tfs.make_new_fig(fnum=1, figsize=fsize)

tick_fs=8
label_fs=10

phs = []
legends = [r'$\mathit{\phi}$ = $\frac{t^2}{2}$',
		   r'$\mathit{\phi}$ = $\frac{t^3}{3}$',
		   r'$\mathit{\phi}$ = $\frac{t^4}{4}$',
		   r'$\mathit{\phi}$ = $e^{-t/(6s)}$',
		   r'$\mathit{\phi}$ = $\mathit{ln}$ t']

ymin = None;		ymax = None
xmin = min(tau);	xmax = max(tau)
for ii in range(0, 5):
	ph, = plot(tau, oe2e[:, ii], '--o', lw=linewidth, ms=markersize, mew=markeredge)
	ph.set_mec(ph.get_color())
	phs.append(ph)

	cymin = min(oe2e[:, ii])
	cymax = max(oe2e[:, ii])

	ymin = cymin if ymin is None else min([ymin, cymin])
	ymax = cymax if ymax is None else max([ymax, cymax])

xran = xmax-xmin;			yran = ymax-ymin
xbuff = 0.01;				ybuff = 0.05;

xlim([xmin-xran*xbuff, xmax+xran*xbuff])
ylim([ymin-yran*ybuff, ymax+yran*ybuff])

xlabel(r'$\tau$ (s)', fontsize=label_fs);
ylabel(r'$\mathit{\phi}$ (rad)', fontsize=label_fs)
title(r'Phase accumulation after 2s')

leg = legend(phs, legends, loc=2, numpoints=1, handlelength=4, handletextpad=0.1)
leg.legendPatch.set_facecolor('none');
tick_params('both', labelsize=tick_fs, top=False, right=False)

tfs.save_figs('../relaxometry/PhaseAccumulationTaus', formats=['eps', 'pdf'])

tight_layout(h_pad=0.35)
subplots_adjust(left=0.09, right=0.98, bottom=0.07, top=0.95)

show();

'''
Second Figure:
Phase accumulation during a CPMG
'''
fnum = tfs.make_new_fig(fnum, figsize=fsize)

ss = arange(0, len(t), 18);

phs = []
legends = []
lstr = r'$\tau$ = {:0.3f} s'
sels = [0, 1, 2, 4, 9]
czorder = len(sels)+1

ymin = None;		ymax = None

for jj, ii in enumerate([0, 1, 2, 4, 9]):
	ph, = plot(t[ss], o[ii, ss], '-');
	
	ph.zorder = czorder
	czorder -= 1

	ph.set_mfc('none')
	ph.set_ms(markersize)
	ph.set_lw(linewidth)

	phs.append(ph)
	legends.append(lstr.format(tau[jj]))

	cymin = min(o[ii, ss])
	cymax = max(o[ii, ss])

	ymin = cymin if ymin is None else min([ymin, cymin])
	ymax = cymax if ymax is None else max([ymax, cymax])


xlabel(r't (s)', fontsize=label_fs);
ylabel(r'$\mathit{\phi}$ (rad)', fontsize=label_fs, labelpad=-1)

yran = ymax-ymin;		ybuff = 0.01
xlim(min(t[ss]), max(t[ss]))
ylim(ymin-ybuff*yran, ymax+ybuff*yran)


title(r'Phase accumulation during a CPMG $\left(\!\ \mathit{\phi} =\ t^{2}\!\ \right)$')

leg = legend(phs, legends, loc=2, handlelength=4)  
leg.legendPatch.set_facecolor('none')

yt = gca().get_yticks();
yt1 = [0]*(len(yt)-1);
for ii in range(0, len(yt)-1):
    yt1[ii] = 0.5*(yt[ii]+yt[ii+1])
gca().set_yticks(yt1);
tick_params('both', labelsize=tick_fs, top=False)

tight_layout(h_pad=0.35)
subplots_adjust(left=0.09, right=0.98, bottom=0.07, top=0.95)

tfs.save_figs('../relaxometry/PhaseAccumulationt', formats=['eps', 'pdf'])

show();
