"""
Plot the results of a quadrature-detected indirect FID.

@author Paul J. Ganssle
"""
import scipy.io;
from scipy import interpolate
from numpy import *;

from matplotlib import rcParams, rc
from matplotlib.pyplot import *;

from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

import matplotlib.gridspec as gridspec

import theme.colors as thc
import theme.figure_settings as tfs

from general_util import *

from analyzedata.fileio.importdata import importdata as idata;
tfs.setup_environment();

fnum = 1

# Read the data into variablesd
out = idata(r'data_sets\quad_fid.mat');
t1 = out['t1'];
cx = mean(out['cx'],1)
cy = mean(out['cy'],1)

t1 = t1[1:];		cx= cx[1:];		cy = cy[1:]
t1 *= 1000					# Convert to ms
cc = cx + 1j*cy;			# Represent as a complex number
cc -= mean(cc)
cc /= max(abs(cc))			# Normalize the data

phase_off = arctan2(imag(cc), real(cc))
phase_off = mean(phase_off[2:10])
cc = -cc									# This is inverted for whatever reason.

shift_phase = False
if shift_phase:
	cc = cc*exp(1j*(phase_off+pi/32+pi))		# Roughly phase shift this bad boy.


# Plot the actual figures
figsize=(6, 3)
dpi = 150
[fig, fnum] = tfs.make_new_fig(fnum=fnum, close_old=True, return_fig=True, 
                               figsize=figsize, dpi=dpi)

xx = linspace(min(t1), max(t1), 1024)
tckr = interpolate.splrep(t1, real(cc), s=0.0075)
tcki = interpolate.splrep(t1, imag(cc), s=0.0075)

yyr = interpolate.splev(xx, tckr)
yyi = interpolate.splev(xx, tcki)

real_color = thc.bc;		imag_color = thc.gc
marker_size=3
lw = 1.25

# Plost basically nothing just to get handles with both things on it.
phr, = plot([0, 0], [0, 1], '-o', color=real_color, mec=real_color, ms=marker_size)
phi, = plot([0, 0], [0, 1], '-o', color=imag_color, mec=imag_color, ms=marker_size)
cla()			# Clear the current axes and start over.

# The real plots
plot(t1, real(cc), 'o', color=real_color, mec=real_color, ms=marker_size)
plot(xx, yyr, '-', lw=lw, color=real_color)

plot(t1, imag(cc), 'o', color=imag_color, mec=imag_color, ms=marker_size)
plot(xx, yyi, '-', lw=lw, color=imag_color)

ymin = min([min(yyr), min(yyi), min(real(cc)), min(imag(cc))]);
ymax = max([max(yyr), max(yyi), max(real(cc)), max(imag(cc))])

yrange = ymax-ymin;			ybuff = 0.05
ylow = ymin-yrange*ybuff;		yhigh = ymin+yrange*ybuff
if ylow < 0:
	ylow = floor(ylow/0.2)*0.2
else:
	ylow = 0

    xlim(ylow, yhigh)

# Labels and annotations and such
legend([phr, phi], ['X Channel', 'Y Channel'], loc='best', 
       fontsize=10, handlelength=3, markerscale=2, numpoints=1)
ylabel('Magnetization (Normalized)')
xlabel('Time (ms)')

gca().set_xticks(arange(0, max(t1)+1, 50))
tick_params(axis='both', labelsize=7, left=False, right=False, top=False)

tight_layout(pad=0.35)

tfs.save_figs('../relaxometry/quadrature_fid_data', formats=['eps', 'pdf'])

show()