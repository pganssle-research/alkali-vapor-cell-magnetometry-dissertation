"""
Generate figures from the decane/water T1-T2 and T2-D spectra.

@author Paul J. Ganssle
"""

from scipy.io import loadmat

import theme.colors as thc
import theme.figure_settings as tfs
import os

from matplotlib.pyplot import *
from numpy import squeeze, array, linspace, transpose
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# General setup
fnum = 0
fsize=(6.5, 6.5)

tfs.setup_environment(figsize=fsize,
                      dpi=150,
                      legend_fs=9,
                      tick_fs=9,
                      label_fs=11,
                      ann_fs=12,
                      n_legend_points=1)

loc = '../relaxometry/'

#
# Load and process data
#

# T1-T2 Data
o = loadmat('data_sets/T1T2DataLite2.mat')
o2 = o['o2']
tau1 = squeeze(o['tau1'])
tau2 = squeeze(o['tau2'])
f = o2['f'][0][0][0][0]
ds = o['ds']
x = squeeze(ds['x'][0][0])
y = squeeze(ds['y'][0][0])
z = squeeze(ds['z'][0][0])

# Get relative sums - this is not information displayed anywhere
av = sum(f, 0)
tv = sum(av)
av /= tv
dsumt1t2 = sum(av[tau1 < 2])
wspan = [x > 2 and x < 6 for x in tau1]
wsumt1t2 = sum(av[array(wspan)])

tau1_t1t2 = tau1;           tau2_t1t2 = tau2
f_t1t2 = 100*f/tv;          tv_t1t2 = tv;

# Diffusion Data
o = loadmat('data_sets/DiffusionT2DataLite2.mat')
o2 = o['o2']
tau1 = squeeze(o['tau1'])
tau2 = squeeze(o['tau2'])
f = o2['f'][0][0][0][0]

av = sum(f, 0)
tv = sum(av)
av /= tv
dsumt2d = sum(av[tau1 < 2])
wspan = [x > 2 and x < 6 for x in tau1]
wsumt2d = sum(av[array(wspan)])

tau1_t2d = tau1;        tau2_t2d = tau2;
f_t2d = 100*f/tv

#
#   Plot the data
#

# T1-T2
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True)

ca = contourf(tau1_t1t2, tau2_t1t2, transpose(f_t1t2), levels=linspace(0.075, 3.5, 250))
xlim(0, 5);
ylim(0, 5);
subplots_adjust(left=0.085, bottom=0.11, right=0.98, top=0.96)

# Labels
xlabel('T$_\mathregular{1}$ (s)')
ylabel('T$_2$')

# Annotations
annote_dict = {'xycoords':'data',
               'textcoords':'data',
               'ha':'center'}

dpos = {'xy':(1.42, 1.36), 'xytext':(1.4, 1.1)}
wpos = {'xy':(3.88, 3.38), 'xytext':(4.1, 3.55)}
use_arrows = False
if use_arrows:
    annote_dict['arrowprops']=dict(arrowstyle='-|>', connectionstyle='arc3,rad=-0.25', fc='k')
    wpos['xytext'] = (2.5,3.8)
    dpos['xytext'] = (2.85,0.75)

wa = gca().annotate('Water', xy=wpos['xy'], xytext=wpos['xytext'], **annote_dict)
da = gca().annotate('Decane', xy=dpos['xy'], xytext=dpos['xytext'], **annote_dict)

ains = inset_axes(gca(), width='5%', height='60%', loc=2)
colorbar(ca, cax=ains, orientation='vertical', 
         ticks=[round(x*10.0)/10.0 for x in linspace(0.075, 3.3, 10)])


tfs.save_figs(os.path.join(loc,'decane-water-t1-t2'), formats=['pdf', 'eps', 'png'])

# T2-D
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True)

ca = contourf(tau1_t2d, tau2_t2d, transpose(f_t2d), linspace(0.05, 2.7, 125), ls=None)
xlim(0, 5);
ylim(0, 5);

subplots_adjust(left=0.085, bottom=0.11, right=0.98, top=0.96)

# Labels
ylabel('$T_{2}$')
xlabel('Self-diffusion Coefficient ($\mathregular{10^{-5} cm^2\!/s}$)')

# Annotations
annote_dict = {'xycoords':'data',
               'textcoords':'data',
               'va':'center'}

dpos = {'xy':(1.9, 1.48), 'xytext':(2.1, 1.375)}
wpos = {'xy':(3.05, 3.45), 'xytext':(2.85, 3.45)}

if use_arrows:
    annote_dict['arrowprops']=dict(arrowstyle='-|>', connectionstyle='arc3,rad=-0.25', fc='k')
    wpos['xytext'] = (1.75, 3.75)
    dpos['xytext'] = (3.35, 0.75)

wa = gca().annotate('Water', xy=wpos['xy'], xytext=wpos['xytext'], ha='right', **annote_dict)
da = gca().annotate('Decane', xy=dpos['xy'], xytext=dpos['xytext'], ha='left', **annote_dict)

ains = inset_axes(gca(), width='5%', height='60%', loc=1)
colorbar(ca, cax=ains, orientation='vertical', 
         ticks=[round(x*10.0)/10.0 for x in linspace(0.1, 4.6, 10)])
ains.yaxis.set_ticks_position('left')

tfs.save_figs(os.path.join(loc,'decane-water-t2-d'), formats=['pdf', 'eps', 'png'])

show()
