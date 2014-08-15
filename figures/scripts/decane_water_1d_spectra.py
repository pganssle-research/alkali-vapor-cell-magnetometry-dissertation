"""
Generate the 1D Laplace spectra for the decane-water data.

@author Paul J. Ganssle
"""

from scipy.io import loadmat
from scipy.interpolate import interp1d

import theme.colors as thc
import theme.figure_settings as tfs
import os

from matplotlib.pyplot import *
from numpy import squeeze, array, linspace, transpose

# General setup
fnum = 0
fsize=(8, 4)

tfs.setup_environment(figsize=fsize,
                      dpi=150,
                      legend_fs=9,
                      tick_fs=8,
                      label_fs=12,
                      ann_fs=12,
                      n_legend_points=1)

loc = '../relaxometry/'

dcol = thc.rc
wcol = thc.bc

# Load the data
o = loadmat('data_sets/DecaneWaterRelaxation.mat')
o1wt1 = o['o1wt1'];     o1wt2 = o['o1wt2'];
o1dt1 = o['o1dt1'];     o1dt2 = o['o1dt2'];
o1dwt1 = o['o1dwt1'];   o1dwt2 = o['o1dwt2'];

twt1 = squeeze(o1wt1['t'][0][0][0]);
tdt1 = squeeze(o1dt1['t'][0][0][0]);
tdwt1 = squeeze(o1dwt1['t'][0][0][0]);

twt2 = squeeze(o1wt2['t'][0][0][0]);
tdt2 = squeeze(o1dt2['t'][0][0][0]);
tdwt2 = squeeze(o1dwt2['t'][0][0][0]);

a_t1 = 4;
a_t2 = 4;

fwt1 = squeeze(o1wt1['f'][0][0][0][a_t1])
fdt1 = squeeze(o1dt1['f'][0][0][0][a_t1])
fdwt1 = squeeze(o1dwt1['f'][0][0][0][a_t1])

fwt2 = squeeze(o1wt2['f'][0][0][0][a_t2-2])
fdt2 = squeeze(o1dt2['f'][0][0][0][a_t2-3])
fdwt2 = squeeze(o1dwt2['f'][0][0][0][a_t2])

# Calculate the area sums
tv1 = sum(fdwt1)
dsum1 = sum(fdwt1[tdwt1 < 2])
wsum1 = sum(fdwt1[tdwt1 >= 2])

tv2 = sum(fdwt2)
dsum2 = sum(fdwt2[tdwt2 < 2])
wsum2 = sum(fdwt2[tdwt2 >= 2])

dfrac1 = dsum1/tv1;        wfrac1 = wsum1/tv1;
dfrac2 = dsum2/tv2;        wfrac2 = wsum2/tv2; 

# Smooth out the data with some spline fitting
tnew = linspace(0.5, 4.5, 1000)
sfdwt1 = interp1d(append(0, tdwt1), append(0, fdwt1), kind='cubic')
sfdt1 = interp1d(append(0, tdt1), append(0, fdt1), kind='cubic')
sfwt1 = interp1d(append(0, twt1), append(0, fwt1), kind='cubic')

sfdwt2 = interp1d(append(0, tdwt2), append(0, fdwt2), kind='cubic')
sfdt2 = interp1d(append(0, tdt2), append(0, fdt2), kind='cubic')
sfwt2 = interp1d(append(0, twt2), append(0, fwt2), kind='cubic')

fdwt1 = sfdwt1(tnew)
fdt1 = sfdt1(tnew)
fwt1 = sfwt1(tnew)

fdwt2 = sfdwt2(tnew)
fdt2 = sfdt2(tnew)
fwt2 = sfwt2(tnew)

tdwt1 = tnew;       tdt1 = tnew;        twt1 = tnew;
tdwt2 = tnew;       tdt2 = tnew;        twt2 = tnew;

# Decane/Water T1
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True)

alpha_level = 0.8;
flw = 1.25;
y_buff = 1.025
fill_between_options = {
  'interpolate':True,
  'alpha':alpha_level,
  'antialiased':True,
  'lw':0
}


subplot(3, 2, 5)     
plot(append(0, tdwt1), append(0, fdwt1), 'k', lw=flw);
fill_between(tdwt1, fdwt1, where=(tdwt1<2.5), facecolor=dcol, **fill_between_options);
fill_between(tdwt1, fdwt1, where=(tdwt1>2.5), facecolor=wcol, **fill_between_options);
xlim(0.1, 4.6)
ylim(0, max(fdwt1)*y_buff)
yticks([])
gca().set_alpha(0)

subplot(3, 2, 1)
plot(append(0, tdt1), append(0, fdt1), 'k', lw=flw);
fill_between(tdt1, fdt1, facecolor=dcol, **fill_between_options);

xlim(0.1, 4.6)
ylim(0, max(fdt1)*y_buff)
yticks([])
gca().set_alpha(0)

subplot(3, 2, 3)
plot(append(0, twt1), append(0, fwt1), 'k', lw=flw);
fill_between(twt1, fwt1, facecolor=wcol, **fill_between_options);

xlim(0.1, 4.6)
ylim(0, max(fwt1)*y_buff)
yticks([])
gca().set_alpha(0)

subplot(3, 2, 5)
xlabel('T$_1$ (s)');

# Decane/Water T2
subplot(3, 2, 6)
plot(append(0, tdwt2), append(0, fdwt2), 'k', lw=flw);
fill_between(tdwt2, fdwt2, where=(tdwt2<2.5), facecolor=dcol, **fill_between_options);
fill_between(tdwt2, fdwt2, where=(tdwt2>2.5), facecolor=wcol, **fill_between_options);

xlim(0.1, 4.6)
ylim(0, max(fdwt2)*y_buff)
yticks([]);
gca().set_alpha(0)

subplot(3, 2, 2)
plot(append(0, tdt2), append(0, fdt2), 'k', lw=flw);
fill_between(tdt2, fdt2, facecolor=dcol, **fill_between_options);
xlim(0.1, 4.6)
ylim(0, max(fdt2)*y_buff)
yticks([]); 
gca().set_alpha(0)

subplot(3, 2, 4)
plot(append(0, twt2), append(0, fwt2), 'k', lw=flw);
xlim(0.1, 4.6)
ylim(0, max(fwt2)*y_buff)
yticks([]); 
fill_between(twt2, fwt2, facecolor=wcol, **fill_between_options);
gca().set_alpha(0)

subplot(3, 2, 6)
xlabel('T$_2$ (s)')

subplots_adjust(left=0.02, bottom=0.125, right=0.98, top=0.98, wspace=0.07, hspace=0.2)


tfs.save_figs(os.path.join(loc,'decane-water-t1-t2-1d'), formats=['pdf', 'eps', 'png'])
show()
