'''
Plot a figure with the CYCLOPS sequence
'''
import scipy.io;
from numpy import *;

from matplotlib import rcParams, rc
from matplotlib.pyplot import *;

import matplotlib.gridspec as gridspec

import theme.colors as thc
import theme.figure_settings as tfs

tfs.setup_environment()

try:
    SaveFigs
except:
    SaveFigs = False

fsize = (7, 4)
(fig, fnum) = tfs.make_new_fig(fnum=1, figsize=fsize, return_fig=True)


## Plot all the spectra
T1 = 1
T2 = 1
freq = 0
np = 2**16

t = linspace(0, 1000, np)
sig = exp(-t/T2)*exp(-1j*2*pi*freq*t)

sr = 1/double(diff(t[0:2]))

f = linspace(-sr/2, sr/2, np)
spec = fftshift((2.0/np)*fft.fft(sig))
spec /= double(max(abs(spec)))    

phase_response = [0, 1, 2];
each_color = [thc.bc, thc.rc, thc.gc]

# Place a gridspec
left_frac = 0.8
bot_frac = 0.70

gsp = gridspec.GridSpec(3, 4,
                       height_ratios=[1, 2, 2],
                       left=0.00,
                       bottom=0.0,
                       right=left_frac,
                       top=bot_frac,
                       wspace=0.02,
                       hspace=0.1)
spec_axes = []
for jj in range(0, 3):
    for ii in range(0, 4):
        ax = subplot(gsp[jj, ii])
        if ii == 0:
            spec_axes.append(ax)

        x_lim = 10
        fsi = find(f >= -x_lim)[0];     fei = find(f <= x_lim)[-1]
        
        cphase = (ii-1)*phase_response[jj]*(pi/2)
        
        
        plot(f[fsi:fei], real(spec[fsi:fei]*exp(1j*cphase)), lw=1.5, color=each_color[jj],
             solid_capstyle='round', solid_joinstyle='round')
        
        xticks([])
        yticks([])
        
        xlim(-x_lim*1.05, x_lim*1.05)
        
        if jj == 0:
            ylim(-0.025, 1.025)
        else:
            ylim(-1.05, 1.05)
        
        ax.axis('off')

## Draw the circles for the receiver/transmitter diagram
# Colors
circ_color = thc.grc
circ_ax_color = thc.grc
dot_color = thc.oc
outline_color='k'

# Line-widths
circ_lw=2.5
circ_ax_lw = 1.5
outline_lw = 1.5

# Geometry
circ_rad = 1.0/4.0
line_rad = circ_rad*1.25
dot_r = circ_rad/5
dot_cr = circ_rad*1.4+dot_r
dot_c_l = (0.5, 0.5+dot_cr)
dot_z = 1

outline_r = dot_r
outline_c_x = [0.5+x for x in [-dot_cr, 0, dot_cr, 0]]
outline_c_y = [0.5+y for y in [0, dot_cr, 0, -dot_cr]]
outline_z = 2

aspect_ratio = (1.0*fsize[1])/(1.0*fsize[0])
row_h = (1.0-bot_frac)*fsize[1]
fill_w = left_frac*fsize[0]
col_w = fill_w/4.0
extra_space = (4.0*(col_w-row_h))/fsize[0]      # In fraction of the figure
lr_marg = (extra_space/8)
wspace = (extra_space/4)

col_f = row_h/fsize[0]
row_f = row_h/fsize[1]

for ii in range(0, 4):
    ax = fig.add_axes([lr_marg+ii*(col_f+wspace), bot_frac, col_f, row_f], frameon=False, axisbg='none')
    
    circ = Circle((0.5,0.5), circ_rad, lw=circ_lw, color=circ_color, fill=False, transform=ax.transAxes)
    dot = Circle(dot_c_l, dot_r, lw=0, fill=True, color=dot_color, zorder=dot_z, transform=ax.transAxes)
    outline = Circle((outline_c_x[ii], outline_c_y[ii]), dot_r, lw=outline_lw, color=outline_color, fill=False, zorder=outline_z, transform=ax.transAxes)


    lx = Line2D([0.5, 0.5], [0.5-line_rad, 0.5+line_rad], lw=circ_ax_lw, color=circ_ax_color, transform=ax.transAxes)
    ly = Line2D([0.5-line_rad, 0.5+line_rad], [0.5, 0.5], lw=circ_ax_lw, color=circ_ax_color, transform=ax.transAxes) 

    ax.add_patch(circ)
    ax.add_patch(dot)
    ax.add_patch(outline)
    ax.add_line(lx)
    ax.add_line(ly)
    ax.axis('off')

## Add the annotations
cell_h = (1-bot_frac)
cell_w = (1-left_frac)
dot_r *= cell_h*1.6

lr_marg = -0.1
tb_marg = 0.25
text_marg = 0.025

rec_line_h = bot_frac+cell_h*(1-tb_marg)-dot_r
trans_line_h = bot_frac+tb_marg*cell_h

ax = fig.add_axes([left_frac+cell_h*lr_marg+dot_r, rec_line_h, dot_r, dot_r*fsize[0]/fsize[1]])
dot = Circle((0.5, 0.5), 0.95/2, lw=0, fill=True, color=dot_color, zorder=dot_z, transform=ax.transAxes)
ax.add_patch(dot)
ax.axis('off')

ax = fig.add_axes([left_frac+cell_h*lr_marg+dot_r, trans_line_h, dot_r, dot_r*fsize[0]/fsize[1]])
outline = Circle((0.5, 0.5), 0.75/2, lw=outline_lw, color=outline_color, fill=False, zorder=outline_z, transform=ax.transAxes)
ax.add_patch(outline)
ax.axis('off')

fig.text(x=left_frac+cell_h*(lr_marg+text_marg)+dot_r*2, y=rec_line_h-0.008 ,
        s='Receiver',
        fontsize=15,
        ha='left',
        va='bottom')

fig.text(x=left_frac+cell_h*(lr_marg+text_marg)+dot_r*2, y=trans_line_h-0.008 ,
        s='Transmitter',
        fontsize=15,
        ha='left',
        va='bottom')

top_frac = (1-bot_frac)
ratios = [0.0, 1.0, 2.0, 2.0]
bottoms = cumsum(ratios)/sum(ratios)
centers = 1-(top_frac+(bot_frac*array([mean([bottoms[ii], bottoms[ii+1]]) for ii in range(0, len(ratios)-1)])))

for ii in range(0, 3):
    # Now the p = annotations.
    fig.text(x=mean([1, left_frac]), y=centers[ii],
             s='$\mathit{{p={:0.0f}}}$'.format(ii),
             fontsize=20,
             ha='center',
             va='center')


tfs.save_figs('../console/CYCLOPS', formats=['eps', 'pdf'])

show()