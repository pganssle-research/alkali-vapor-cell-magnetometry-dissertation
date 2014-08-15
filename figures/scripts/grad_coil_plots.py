"""
Generate plots for the gradient coil calibration

@author Paul J. Ganssle
"""
from __future__ import division

from scipy.io import loadmat
from numpy import *;
from os.path import join as pathjoin

from matplotlib import rcParams, rc
from matplotlib.pyplot import *;

from matplotlib.patches import Rectangle, FancyArrowPatch, Arc
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D


from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import theme.colors as thc
import theme.figure_settings as tfs

from matplotlib.colors import ColorConverter

import argparse

from general_util import maxf, minf, find_where

parser = argparse.ArgumentParser(description="Generate pulse sequence acquisitions.")
parser.add_argument('-s', '--save', help='Save the figures', action='store_true')
parser.add_argument('-c', '--close', help='Close the figures when done', action='store_true')
args = parser.parse_args()

SaveFigs = args.save
CloseFigs = args.close

tfs.setup_environment()

loc = '../coils'
make_path = lambda name: pathjoin(loc, name)
formats = ['png', 'eps']

figsize=(13, 3)
dpi = 150

# Load data
fname = 'data_sets/2012-08-16-GradCoilImage.mat'
o = loadmat(fname)

# Gradient coil fields
xx = squeeze(o['xpfz'])*25.4            # Conver to mm
zz = squeeze(o['zpfz'])*25.4
Gzf = squeeze(o['Gzpfz'])

# Gradient coil FID, simulations
gpf1r = squeeze(o['gpf1r'])
fn = squeeze(o['fn'])
G_str = squeeze(o['G_str'])
Z_str = squeeze(o['z_str'])

# Gradient coil, image.
f = squeeze(o['f2'])
sg = o['s3']
sgp = sg[:, 0]
sgn = sg[:, 1]

# Calibrations
z_cal_s = Z_str;                 # G
z_cal_e = 0.5075;               # G
gammah = 4257;                  # Hz/G
res_freq_s = gammah*z_cal_s;    # Hz
res_freq_e = gammah*z_cal_e;    # Hz

f -= res_freq_e

# Plot gradient responses.
[fig, fnum] = tfs.make_new_fig(fnum=1, close_old=True, return_fig=True,
                             figsize=figsize, dpi=dpi)

# Find the distance from the center frequency to each of the tops.
w1 = 2;         w2 = 1;
ncol= w1+w2

ax1 = subplot2grid((1, ncol), (0,0), colspan=w1)
iip = argmax(abs(diff(sgp)))+1
iin = argmax(abs(diff(sgn)))-1

ymax = max([max(sgp), max(sgn)])
ybuff = 0.2

yr = [0, ymax*(1+ybuff)]
dash_opts = {
    'linestyle':'--',
    'color':'k',
    'lw':1
}

# Plot the data
pp, = plot(f, sgp, color=thc.bc)
pn, = plot(f, sgn, color=thc.rc)

# Plot dashed lines at the inflection points
plot([f[iin]]*2, yr, **dash_opts)
plot([f[iip]]*2, yr, **dash_opts)
plot([0]*2, yr, **dash_opts)

# Arrows with text
arrow_kwargs = {
    'arrowstyle':'<|-|>',
    'fc':'k',
    'mutation_scale':25,
    'transform':gca().transData,
    'lw':1.25
}

text_kwargs = {'ha':'center',
    'fontsize':8,
    'transform':gca().transData
}

ah = ymax*(1+ybuff/2)
arrow_p = FancyArrowPatch((0, ah), (f[iip], ah), **arrow_kwargs)
arrow_n = FancyArrowPatch((0, ah), (f[iin], ah), **arrow_kwargs)


text(mean([0, f[iip]]), ah*(1+0.005), '{:0.3f} mm'.format(5.066),
     va='bottom',
     **text_kwargs)
text(mean([0, f[iip]]), ah*(1-0.0175), '{:0.1f} Hz'.format(f[iip]),
     va='top',
     **text_kwargs)

text(mean([0, f[iin]]), ah*(1+0.005), '{:0.3f} mm'.format(5.066),
     va='bottom',
     **text_kwargs)
text(mean([0, f[iin]]), ah*(1-0.0175), '{:0.1f} Hz'.format(f[iin]),
     va='top',
     **text_kwargs)

gca().add_patch(arrow_p)
gca().add_patch(arrow_n)

xlim([-1500, 1500])
ylim(yr)

tick_params(top=False, left=False, right=False)

ylabel('Magnetization (pT)')
xlabel('Frequency Offset (Hz)')

IG = 7.5*0.0658;  
legend([pp, pn], 
       ['I$_G$ = + {:0.2f} A'.format(IG), 'I$_G$ = {:0.2f} A'.format(-IG)], 
       fontsize=11,
       loc='best')

text(0.01, 0.975, 'a.)', 
        transform=ax1.transAxes,
        ha='left',
        va='top',
        fontsize=13)

# Now compare the positive data to simulations
ax2 = subplot2grid((1, ncol), (0,w1), colspan=w2)
gpf1r = fftshift(abs(gpf1r))
gpf1r *= maxf(sgp)/maxf(gpf1r)

fn /= gammah*G_str/10
fn += 3.24
f *= 10*0.50658/f[iip]

ps, = plot(fn, gpf1r, '--', color=thc.gc)
pe, = plot(f, sgp, color=thc.bc)

xlim([-15, 15])

legend([ps, pe], ['Simulation', 'Experiment'], loc='best')

text(0.975, 0.975, 'b.)', 
        transform=ax2.transAxes,
        ha='right',
        va='top',
        fontsize=13)

tick_params(left=False, right=False, top=False)

xlabel('Longitudinal Displacement (mm)')

tight_layout()

if SaveFigs:
    tfs.save_figs(make_path('GradientImage'), formats=formats)

# Plot the gradient strength along the x axis and y axis
tfs.setup_environment(tick_fs=8)

[fig, fnum] = tfs.make_new_fig(fnum=fnum, close_old=True, return_fig=True,
                             figsize=(13, 7), dpi=dpi)

r0 = 2.1
t_bot = -25.4*(0.125-0.0125)-0.38

sa = find_where(lambda x: x >= t_bot and x <= t_bot+7.9, zz)
sa.insert(0, sa[0]-1)                               # For the diff.

ra = find_where(lambda x: x >= -2.5 and x <= 2.5, xx)

Gzsa = Gzf[:, sa]                   
Gzsd = diff(Gzsa)/diff(zz[sa]/10)
zd = zz[sa[1:]]
xr = xx[ra]
Gzsdr = Gzsd[ra, :]

ca = contourf(xr, zd, Gzsdr.T, 250)

ax = gca()

xlabel('Radial Distance from Coil Center (mm)')
ylabel('Longitudinal Distance from Coil Center (mm)')
ylim(min(zd), max(zd))

tick_params(top=False, bottom=False, left=False, right=False)

tight_layout()

'''
linel = Line2D([-r0, -r0],
               [max(zd), t_bot+r0], transform=ax.transData, color='k')

liner = Line2D([r0, r0],
               [max(zd), t_bot+r0], transform=ax.transData, color='k')

barc = Arc((0, t_bot+r0), height=2*r0, width=2*r0, theta1=180, theta2=0, transform=ax.transData)

ax.add_line(linel)
ax.add_line(liner)
ax.add_patch(barc)
'''

ax_z = inset_axes(ax, width=2.5, height=2, bbox_to_anchor=(0.045, 0.975), 
                  bbox_transform=ax.transAxes, loc=2)
cc = ColorConverter()
ax_z.set_axis_bgcolor(cc.to_rgba('w', alpha=0.6))
xz = argmin(abs(xx))

plot(zd, Gzsd[xz, :], color='k', lw=1.25)

xlim(min(zd), max(zd))
ylim(0.6, 0.65)
yticks(arange(0.6, 0.65, 0.01))
ymin = min(Gzsd[xz, :])
ymax = max(Gzsd[xz, :])

text(0.98, 0.97, 'r = 0 mm', fontsize=8, fontweight='bold', 
     transform=ax_z.transAxes,
     ha='right',
     va='top')

text(-0.125, 0.5, 'Coil Response (G/cm$\cdot$A)', fontsize=8, 
     transform=ax_z.transAxes,
     ha='center',
     va='center',
     rotation='vertical')
xlabel('z (mm)', fontsize=7)

tick_params(top=False, bottom=False, left=False, right=False, labelsize=7, labelcolor='k')

ax_r = inset_axes(ax, width=2.5, height=1.5, bbox_to_anchor=(0.775, 0.325), 
                  bbox_transform=ax.transAxes, loc=2)

ax_r.set_axis_bgcolor(cc.to_rgba('w', alpha=0.6))
zz0 = argmin(abs(zd))

plot(xx, Gzsd[:, zz0], color='k', lw=1.25)

xlim(min(xr), max(xr))
ylim(0.647, 0.653)
#yticks(arange(0.64, 0.67, 0.01))
ymin = min(Gzsd[:, zz0])
ymax = max(Gzsd[:, zz0])

text(0.98, 0.97, 'z = 0 mm', fontsize=8, fontweight='bold', 
     transform=ax_r.transAxes,
     ha='right',
     va='top')

text(-0.125, 0.5, 'Coil Response (G/cm$\cdot$A)', fontsize=8, 
     transform=ax_r.transAxes,
     ha='center',
     va='center',
     rotation='vertical')
xlabel('r (mm)', fontsize=7)

#ax_r.xaxis.set_ticks_position('top')
#ax_r.xaxis.set_label_position('top')

tick_params(top=False, bottom=False, left=False, right=False, labelsize=6, labelcolor='k')

ains = inset_axes(ax, width='2%', height='40%', loc=1)

ming = minf(Gzsdr)
maxg = maxf(Gzsdr)

mingr = round(ming*100)/100
maxgr = round(maxg*100)/100

cb = colorbar(ca, cax=ains, orientation='vertical', 
         ticks=arange(mingr, maxgr, 0.005))

cb.set_label('Coil Response (G/cm$\cdot$A)')

ains.yaxis.set_ticks_position('left')
ains.yaxis.set_label_position('left')

if SaveFigs:
    tfs.save_figs(make_path('GradientMap'), formats=formats)

show()

if CloseFigs:
    tfs.close_all_figs();