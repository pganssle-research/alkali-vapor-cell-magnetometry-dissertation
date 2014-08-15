"""
Generate pulse sequence diagrams for demonstration purposes.

@author Paul J. Ganssle
"""
import scipy.io;
from numpy import *;

from matplotlib import rcParams, rc
from matplotlib.pyplot import *;

from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D

import matplotlib.gridspec as gridspec

import theme.colors as thc
import theme.figure_settings as tfs

import argparse

parser = argparse.ArgumentParser(description="Generate pulse sequence acquisitions.")
parser.add_argument('-s', '--save', help='Save the figures', action='store_true')
args = parser.parse_args()

SaveFigs = args.save

tfs.setup_environment()

figsize=(13, 5)
dpi = 150
[fig, fnum] = tfs.make_new_fig(fnum=1, close_old=True, return_fig=True,
                             figsize=figsize, dpi=dpi)
ax = gca()

subplots_adjust(left=0.025, right=0.975, top=0.99, bottom=0.075)

"""
Functions for creating figures
"""
T2=3.0 
T2s=0.0625
t_max=6.25
tau=0.6
freq=10.0
exc_phase=0.0
ref_phase=90.0

pulse_width=0.02
pulse_margin=0.025
pulse_fc=thc.bc
pulse_ec='k'

echo_color=thc.rc
echo_lw=1.5

pulse_lw=1.0
base_lw=2

fill_alpha=0.5

margins=[None, None, None, None]
margin_bottom=0.05
margin_top=0.15
margin_left=None
margin_right=None

baseline_height=None

echo_h_frac=1.4

annotate_pulses=True
ShowPulses=True
ShowBaseline=True
ShowEchoes=True
AnnotatePulses=True

npulses=5

"""
Single spin echo sequence with constant time.
"""
fig = ax.get_figure()
trans = ax.transAxes
ax.axis('off')

exc_phase *= pi/180.0                           # Degrees to radians
ref_phase *= pi/180.0

# Margins can be specified either with a list [bot, top, left, right] or 
# with specific keyword arguments. The keyword arguments override the list items.
if len(margins) < 4:
    for ii in range(len(margins), 4):
        margins.append(None)

default_margins = [0.0, 0.0, 0.0, 0.0]
for ii, margin in enumerate([margin_bottom, margin_top, margin_left, margin_right]):
    if margin is not None:
        margins[ii] = margin
    elif margins[ii] is None:
        margins[ii] = default_margins[ii]

(margin_bottom, margin_top, margin_left, margin_right) = margins


# Setup a spectrum for the initial FID and the echoes
height_area = 1-margin_top-margin_bottom        # These two define the working area                     
timespan_width = (1-2*pulse_margin-margin_left-margin_right)/(1-pulse_width/4)

pulse_wr = (pulse_width)*timespan_width         # Refocusing pulse width.
pulse_we = pulse_wr/2                           # Excitation pulse width

pulse_te = pulse_we*t_max                       # Pulse duration (time)     
pulse_tr = pulse_wr*t_max

np_f = 2**14
np_e = 2**15+1

t_f = linspace(pulse_te/2.0, tau-pulse_te/2.0, np_f)
t_f = linspace(0, tau, np_f)
t_e = linspace(-tau+pulse_tr/2, t_max-2*tau, np_e)
t_e_o = t_e

decay_curve = lambda t, ph: exp(-abs(t)/T2s)*exp(-1j*(2*pi*freq*t+ph))
f_f = decay_curve(t_f, exc_phase)
f_e = decay_curve(t_e, exc_phase)

# Adjust the phases here.
#phase_offset = arctan2(imag(f_f[0]), real(f_f[0]))
phase_offset = pi
f_e *= exp(1j*(2*(ref_phase-exc_phase)-phase_offset))
#f_f *= exp(-1j*(phase_offset))

f_e = real(f_e)
f_f = real(f_f)

threshhold = 0.0015
lif = None
for ii in reversed(range(0, len(f_f))):     # Find the last index above the threshhold
    if abs(f_f[ii]) > threshhold:
        lif = ii
        break

fi = 0;     lie = -1
for ii in range(0, len(f_e)):           # Find the first AND last index above the threshhold
    if abs(f_e[ii]) > threshhold:
        fie = ii
        break

for ii in reversed(range(0, len(f_e))):
    if abs(f_e[ii]) > threshhold:
        lie = ii
        break

## Calculate the geometry involved, all in transAxes units.
span_top = 1-margin_top
span_bottom = margin_bottom
span_left = margin_left
span_right = 1-margin_right

fid_max = max([max(f_e), max(f_f)])
fid_min = min([0, min(f_e), min(f_f)])          # Fix the absolute minimum at zero.
fid_range = fid_max-fid_min                     

if baseline_height is None:
    if fid_min < 0:                                 # This is the height of the base line.
        baseline_height = span_top/(1+(fid_range/(echo_h_frac*abs(fid_min))))
    else:
        baseline_height = 0

    baseline_height *= height_area
    baseline_height += span_bottom

pulse_height = span_top-baseline_height
echo_height = pulse_height*echo_h_frac

tau_width = (tau/t_max)*timespan_width          # Convert tau to width in axis units.

fid_loc = pulse_margin+pulse_we                 # Where the excitation pulse ends
exc_loc = fid_loc-pulse_we                      # Where the excitation pulse starts

ref_loc = exc_loc+tau_width+pulse_we/2-pulse_wr/2.0     # Where the refocusing pulse starts
echo_loc = ref_loc+pulse_wr                     # Where the echo starts.

ip_span = ref_loc-fid_loc                       # Span between the excitation and refocus pulse
ap_span = timespan_width-echo_loc               # Span the remainder of the figure
                                    
## Create the shapes
# Keep track of all these zorders. Line is on the bottom, then FIDs, then pulses on top
line_zorder = 1
fid_zorder = 2
pulse_zorder = 3

if ShowBaseline:
    # Generate the pulses and baseline first.
    baseline = Line2D([span_left, span_right],
                      [baseline_height]*2,
                      lw=base_lw,
                      color='k',
                      zorder=line_zorder,
                      clip_on=False,
                      transform=trans)

if ShowPulses:
    exc_pulse = Rectangle((exc_loc, baseline_height),
                          pulse_we,
                          pulse_height,
                          fc=pulse_fc,
                          ec=pulse_ec,
                          lw=pulse_lw,
                          zorder=pulse_zorder,
                          clip_on=False,
                          transform=trans)

    for ii in range(0, npulses):
        ref_pulse = Rectangle((ref_loc+2*ii*tau_width, baseline_height),
                              pulse_wr,
                              pulse_height,
                              fc=pulse_fc,
                              ec=pulse_ec,
                              lw=pulse_lw,
                              zorder=pulse_zorder,
                              clip_on=False,
                              transform=trans)

        ax.add_patch(ref_pulse)

    # Now add them to the axis.
    ax.add_patch(exc_pulse)
    ax.add_line(baseline)

## Now add in the echo decays
# Convert the data into transAxes and scale them right.
if ShowEchoes:
    f_f /= fid_range/echo_height                    # Normalize to echo height.
    f_f += baseline_height                          # Offset so the zero point is at the baseline

    f_e_max = copy(f_e)
    
    f_e /= fid_range/echo_height
    f_e += baseline_height

    t_f = linspace(fid_loc, ref_loc, len(f_f))
    t_e_span = (max(t_e)-min(t_e))*timespan_width/t_max
    t_e = linspace(echo_loc, echo_loc+t_e_span, len(f_e))
    
    # Plot the FIDs now
    echo_ls = '-'
    echo_line = {
                 'lw':echo_lw,
                 'zorder':fid_zorder,
                 'color':echo_color,
                 'transform':trans,
                 'clip_on':False
                }

    echo_fill = {
                 'lw' : 0,
                 'zorder':fid_zorder-1,
                 'color':echo_color,
                 'transform':trans,
                 'clip_on':False,
                 'alpha':fill_alpha
                 }

    # Remove the zeros
    threshhold = 0.005

    if lif is not None:
        t_f = t_f[0:lif]
        f_f = f_f[0:lif]

    t_e = t_e[fie:lie]
    f_e = f_e[fie:lie]

    plot(t_f, f_f, echo_ls, **echo_line)
    fill_between(t_f, baseline_height, f_f, **echo_fill)

    f_e_max = f_e_max[fie:lie]
    t_e_o = t_e_o[fie:lie]
    
    tcs = []
    adjs = []

    for ii in range(0, npulses):
        f_e = copy(f_e_max)*exp((-2*(ii+1)*tau)/T2)

        f_e /= fid_range/echo_height
        f_e += baseline_height

        am = argmax(f_e)
        tcs.append(t_e[am]+2*ii*tau_width)
        adjs.append(f_e[am])


        plot(t_e+2*ii*tau_width, f_e, echo_ls, **echo_line)
        fill_between(t_e+2*ii*tau_width, baseline_height, f_e, **echo_fill)

    plot(tcs, adjs, 'o', color='k', ms=5, zorder=4)

    t = linspace(0, 2*tau*npulses, 200)
    m = exp(-t/T2)
    
    m /= fid_range/echo_height
    m += baseline_height    
    t = linspace(fid_loc, fid_loc+2*npulses*tau_width, 200)
    plot(t, m, '--', color='k', lw=2.0, zorder=3)
    xlim([0, 1])
    ylim([0, 1])

    ## Annotate the pulses
    if AnnotatePulses:
        pass

# Draw lines below it
ta2_ff = 0.00092235*0.375           # Fudge factor


ta1_s = exc_loc+pulse_we/2;     ta1_e = ta1_s + tau_width
ta2_s = ref_loc+pulse_wr/2;     ta2_e = ta2_s + tau_width + ta2_ff
ta3_s = ta2_e;                  ta3_e = ta3_s + 2*tau_width 

ah = baseline_height-0.15
arrow_ms = 50;          arrow_lw = 2.0
arrow_kwargs = {'arrowstyle':'<|-|>',
                'fc':'k',
                'mutation_scale':arrow_ms,
                'transform':ax.transAxes,
                'zorder':4,
                'lw':arrow_lw}
tau_arrow1 = FancyArrowPatch((ta1_s, ah), (ta1_e, ah), **arrow_kwargs)
tau_arrow2 = FancyArrowPatch((ta2_s, ah), (ta2_e, ah), **arrow_kwargs)
tau_arrow3 = FancyArrowPatch((ta3_s, ah), (ta3_e, ah), **arrow_kwargs)

dash_y = [ah, pulse_height+baseline_height]
dash_kwargs = {
    'lw':1.2,
    'color':'k',
    'clip_on':False,
    'dashes':[5,7.5,5,7.5],
    'transform':ax.transAxes
}
dashs = Line2D([ta1_s, ta1_s], dash_y, zorder=4, **dash_kwargs)
dashm = Line2D([ta2_s, ta2_s], dash_y, zorder=4, **dash_kwargs)
dashe = Line2D([ta2_e, ta2_e], dash_y, zorder=-0, **dash_kwargs)
dashe2 = Line2D([ta3_e, ta3_e], dash_y, zorder=0, **dash_kwargs)

ax.add_patch(tau_arrow1)
ax.add_patch(tau_arrow2)
ax.add_patch(tau_arrow3)
ax.add_line(dashe2)
ax.add_line(dashs)
ax.add_line(dashm)
ax.add_line(dashe)

# Annotate the arrows
txt1x = mean([ta1_s, ta1_e])
txt2x = mean([ta2_s, ta2_e])
txt3x = mean([ta3_s, ta3_e])

txth = ah-0.05
ptxth = baseline_height+pulse_height+0.075

# Fudge factors because we're adding subscripts
p2ff = 0.006
p2hff = 0.002
pff = 0.005
phff = 0.015

txtkwargs = {
    'fontsize':20,
    'va':'center',
    'ha':'center',
    'transform':ax.transAxes
}

#rcParams['text.usetex']=True
text(txt1x, txth, r'$\tau$',**txtkwargs)
text(txt2x, txth, r'$\tau$', **txtkwargs)
text(txt3x, txth, r'2$\tau$', **txtkwargs)

text(exc_loc+pulse_we/2+p2ff, ptxth-p2hff, r'$\frac{\pi}{2}_{x}$', **txtkwargs)

for ii in range(0, npulses):
    text(ref_loc+pulse_wr/2+pff+2*ii*tau_width, ptxth-phff, r'$\pi_y$', **txtkwargs)

# These are for
#text(exc_loc+pulse_we/2, ptxth, r'$\frac{\pi}{2}$', color='r', **txtkwargs)
#text(ref_loc+pulse_wr/2, ptxth, r'$\pi$', color='r', **txtkwargs)

if SaveFigs:
    tfs.save_figs('../relaxometry/cpmg_sequence_demo',
                  formats=['eps', 'png', 'pdf'])

show()


