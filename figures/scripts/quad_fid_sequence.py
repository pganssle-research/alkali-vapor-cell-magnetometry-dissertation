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

qchans = [0]*4
qacq_chan = 3;  qby_chan = 2;   qbx_chan = 1;   qbz_chan = 0;
qchans[qbx_chan] = bx;
qchans[qby_chan] = by;
qchans[qbz_chan] = bz;
qchans[qacq_chan] = ac;

# Acquisition pulse sequence
qfid = PulseSequence(qchans)                          # Quadrature FID acquisition

# Have to build the sequence in order. Start with a margin
lmarg_dur = 20; rmarg_dur = 30
pi_dur = 30
rd_dur = 10
acq_dur = 100
shuttle_dur=400
tau_dur=100

# The labels
rcParams['text.usetex']=False
amsmath = r'\usepackage{xfrac}'
if amsmath not in rcParams['text.latex.preamble']:
    rcParams['text.latex.preamble'].append(amsmath)
label90 = r'$\sfrac{\pi}{2}$' if rcParams['text.usetex'] else r'$\frac{\pi}{2}$'
label180 = r'$\pi$'
label_tshut = r'$t_{shuttle}$'
label_tau = r'$\tau$'
label_acq = r'$t_m$'

### Indirect FID acquisition - quadrature mode
# The pulses
off_stage = qfid.make_stage(None)
x_stage = qfid.make_stage([qbx_chan])
y_stage = qfid.make_stage([qby_chan])
z_stage = qfid.make_stage([qbz_chan])
acq_stage = qfid.make_stage([qacq_chan])

qfid.add_stage(duration=lmarg_dur, states=off_stage)
qfid.add_stage(duration=pi_dur/2, states=x_stage,
            label=label90, labelchan=qbx_chan)
qfid.add_stage(duration=tau_dur, states=z_stage,
            label=label_tau, labelchan=qbz_chan)
qfid.add_stage(duration=pi_dur/2, states=x_stage,
             label=label90, labelchan=qbx_chan)
sec_start = qfid.add_stage(duration=rd_dur, states=off_stage)
qfid.add_stage(duration=acq_dur, states=acq_stage,
               label=label_acq, labelchan=qacq_chan)
qfid.add_stage(duration=rd_dur, states=off_stage)
qfid.add_stage(duration=pi_dur, states=y_stage,
            label=label180, labelchan=qby_chan)
qfid.add_stage(duration=rd_dur, states=off_stage)
qfid.add_stage(duration=acq_dur, states=acq_stage,
               label=label_acq, labelchan=qacq_chan)
qfid.add_stage(duration=rd_dur, states=off_stage)
qfid.add_stage(duration=pi_dur/2, states=x_stage,
            label=label90, labelchan=qbx_chan)
sec_end = qfid.add_stage(duration=rd_dur, states=off_stage)
qfid.add_stage(duration=rmarg_dur, states=off_stage)

# Repeated sections
rep_sec = qfid.add_section(sec_start, label='N', color=thc.rc, lw=1.5,
                        offset_above=0.01, offset_below=0.0175,
                        offset_before=0.25, offset_after=0.35,
                        text_offset_below=0.0175, text_offset_right=0.015)
qfid.end_section(sec_end, section=rep_sec)

qfid.add_line(1, label='Sample Drop', color=thc.rc, lw=1.75, 
            offset_before=0.005, offset_below=0.00, 
            offset_above=0.01, text_offset_below=0.02,
            label_fs=8)

figsize=(6.5, 3)
dpi = 150
[fig, fnum] = tfs.make_new_fig(fnum=fnum, close_old=True, return_fig=True,
                             figsize=figsize, dpi=dpi)


make_pulse_sequence(fig, qfid, 
                    chan_label_weight=0.04, chan_label_right=False)

if SaveFigs:
    tfs.save_figs('../relaxometry/quad_fid_acq_sequence_diagram',
                  formats=['eps', 'pdf', 'png'])

show()