"""
Generate a pulse sequence for acquisitions.
"""

import sequence_diagrams
try:
    SD_loaded

    reload(sequence_diagrams)
except:
    SD_loaded = True

from sequence_diagrams import *
import theme.colors as thc
import theme.figure_settings as tfs

import argparse

parser = argparse.ArgumentParser(description="Generate pulse sequence acquisitions.")
parser.add_argument('-s', '--save', help='Save the figures', action='store_true')
args = parser.parse_args()

SaveFigs = args.save

tfs.setup_environment()

lweight = 0.1
bx = PulseChannel(label=r'B$_\mathsf{X}$', color=thc.gc, labelweight=lweight,
                  label_pos='center', label_fs=12, clabel_fs=15,
                  lw=1)
bz = PulseChannel(label=r'B$_\mathsf{Z}$', color=thc.bc, labelweight=lweight,
                  label_pos='center', label_fs=16, clabel_fs=15,
                  lw=1)
ac = PulseChannel(label=r'Acq.', color=thc.dgrc, labelweight=lweight,
                  label_pos='center', label_fs=12, clabel_fs=15,
                  lw=1)

chans = [0]*3
acq_chan = 2;  bx_chan = 1;     bz_chan = 0;
chans[bx_chan] = bx;
chans[bz_chan] = bz;
chans[acq_chan] = ac;

# Pulse Sequence
t1 = PulseSequence(chans)

# Have to build the sequence in order. Start with a margin
lmarg_dur = 20; rmarg_dur = 30
pi_dur = 25
rd_dur = 10
acq_dur = 100
tau_dur=300

# The labels
rcParams['text.usetex']=False
amsmath = r'\usepackage{xfrac}'
if amsmath not in rcParams['text.latex.preamble']:
    rcParams['text.latex.preamble'].append(amsmath)
label90 = r'$\sfrac{\pi}{2}$' if rcParams['text.usetex'] else r'$\frac{\pi}{2}$'
label180 = r'$\pi$'
label_tshut = r'$t_{shuttle}$'
label_tau = r'$\tau$'
labelacq = r'$t_{m}$'

# The pulses
off_stage = t1.make_stage(None)
x_stage = t1.make_stage([bx_chan])
z_stage = t1.make_stage([bz_chan])
acq_stage = t1.make_stage([acq_chan])

t1.add_stage(duration=lmarg_dur, states=off_stage)
t1.add_stage(duration=tau_dur, states=z_stage,
            label=label_tau, labelchan=bz_chan)
t1.add_stage(duration=pi_dur, states=x_stage,
            label=label180, labelchan=bx_chan)
sec_start = t1.add_stage(duration=rd_dur, states=off_stage)
t1.add_stage(duration=acq_dur, states=acq_stage,
             label=labelacq, labelchan=acq_chan)
t1.add_stage(duration=rd_dur, states=off_stage)
t1.add_stage(duration=pi_dur, states=x_stage,
            label=label180, labelchan=bx_chan)
t1.add_stage(duration=rd_dur, states=off_stage)
t1.add_stage(duration=acq_dur, states=acq_stage,
             label=labelacq, labelchan=acq_chan)
sec_end = t1.add_stage(duration=rd_dur, states=off_stage)
t1.add_stage(duration=rmarg_dur, states=off_stage)

# Repeated sections
rep_sec = t1.add_section(sec_start,
                         label='N',
                         color=thc.rc,
                         lw=1.5,
                         offset_before=0.5,
                         offset_after=0.5,
                         offset_below=0.0175,
                         offset_above=0.01,
                         lip_frac=0.04,
                         text_offset_below=0.025,
                         text_offset_right=0.015)
t1.end_section(sec_end, section=rep_sec)

t1.add_line(1, label='Sample Drop', color=thc.rc, lw=1.75, 
            offset_before=0.01, offset_below=0.00, 
            offset_above=0.01, text_offset_below=0.025,
            label_fs=8)

figsize=(7, 2.5)
dpi = 150
[fig, fnum] = tfs.make_new_fig(fnum=1, close_old=True, return_fig=True,
                             figsize=figsize, dpi=dpi)

make_pulse_sequence(fig, t1,
                    chan_label_weight=0.05,
                    chan_label_right=False)

if SaveFigs:
    tfs.save_figs('../relaxometry/t1_sequence_diagram',
                  formats=['eps', 'png', 'pdf'])

show()