"""
Generate a pulse sequence for FID acquisition.
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
fnum = 0

standard_options = {'labelweight':0.01,
                    'label_pos':'center',
                    'clabel_fs':16,
                    'lw':1.5}
bx = PulseChannel(label=r'B$_\mathsf{X}$', color=thc.gc, label_fs=12, **standard_options)
by = PulseChannel(label=r'B$_\mathsf{Y}$', color=thc.oc, label_fs=16, **standard_options)
bz = PulseChannel(label=r'B$_\mathsf{Z}$', color=thc.bc, label_fs=16, **standard_options)
gz = PulseChannel(label=r'G$_\mathsf{Z}$', color=thc.pc, label_fs=14, **standard_options)
ac = PulseChannel(label=r'Acq.', color=thc.dgrc, label_fs=14, **standard_options)

chans = [0]*3
acq_chan = 2;  bx_chan = 1;     bz_chan = 0;
chans[bx_chan] = bx;
chans[bz_chan] = bz;
chans[acq_chan] = ac;

# Pulse Sequence
fid = PulseSequence(chans)                            # FID Acqusition

# Have to build the sequence in order. Start with a margin
lmarg_dur = 20; rmarg_dur = 30
pi_dur = 30
rd_dur = 10
acq_dur = 100
tau_dur=100
shuttle_dur=150

# The labels
rcParams['text.usetex']=False
amsmath = r'\usepackage{xfrac}'
if amsmath not in rcParams['text.latex.preamble']:
    rcParams['text.latex.preamble'].append(amsmath)
label90 = r'$\sfrac{\pi}{2}$' if rcParams['text.usetex'] else r'$\frac{\pi}{2}$'
label180 = r'$\pi$'
label_tshut = r'$t_{shuttle}$'
label_tau = r'$\tau$'
label_acq = r'$t_{m}$'

#### Indirect FID acquisition.

# The pulses
off_stage = fid.make_stage(None)
x_stage = fid.make_stage([bx_chan])
z_stage = fid.make_stage([bz_chan])
acq_stage = fid.make_stage([acq_chan])

fid.add_stage(duration=lmarg_dur, states=off_stage)
#fid.add_stage(duration=shuttle_dur, states=z_stage,
#            label=label_tshut, labelchan=bz_chan, label_fs=14)
fid.add_stage(duration=pi_dur/2, states=x_stage,
            label=label90, labelchan=bx_chan)
fid.add_stage(duration=tau_dur, states=z_stage,
            label=label_tau, labelchan=bz_chan)
fid.add_stage(duration=pi_dur/2, states=x_stage,
            label=label90, labelchan=bx_chan)
sec_start = fid.add_stage(duration=rd_dur, states=off_stage)
fid.add_stage(duration=acq_dur, states=acq_stage,
              label=label_acq, labelchan=acq_chan)
fid.add_stage(duration=rd_dur, states=off_stage)
fid.add_stage(duration=pi_dur, states=x_stage,
            label=label180, labelchan=bx_chan)
fid.add_stage(duration=rd_dur, states=off_stage)
fid.add_stage(duration=acq_dur, states=acq_stage,
              label=label_acq, labelchan=acq_chan)
sec_end = fid.add_stage(duration=rd_dur, states=off_stage)
fid.add_stage(duration=rmarg_dur, states=off_stage)

# Repeated sections
rep_sec = fid.add_section(sec_start, 
                        label='N', 
                        color=thc.rc, 
                        lw=1.5,
                        offset_before=0.35, 
                        offset_after=0.35,
                        text_offset_below=0.0175, 
                        text_offset_right=0.015)
fid.end_section(sec_end, section=rep_sec)

fid.add_line(1, 
            label='Sample Drop', 
            color=thc.rc, 
            lw=1.75, 
            offset_before=0.005, 
            offset_below=0.00, 
            offset_above=0.01, text_offset_below=0.025,
            label_fs=8)

figsize=(7, 2.5)
dpi = 150
[fig, fnum] = tfs.make_new_fig(fnum=fnum, close_old=True, return_fig=True,
                             figsize=figsize, dpi=dpi)

make_pulse_sequence(fig, fid,                    
                    chan_label_weight=0.05,
                    chan_label_right=False)
if SaveFigs:
    tfs.save_figs('../relaxometry/fid_acq_sequence_diagram',
                  formats=['eps', 'pdf', 'png'])

show()