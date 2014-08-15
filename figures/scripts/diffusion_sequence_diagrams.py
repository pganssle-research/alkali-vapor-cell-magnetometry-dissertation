"""
Generate two diffusion figures - two-pulse and 4-pulse.
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
gz = PulseChannel(label=r'G$_\mathsf{Z}$', color=thc.pc, label_fs=14, **standard_options)
ac = PulseChannel(label=r'Acq.', color=thc.dgrc, label_fs=14, **standard_options)
standard_options['label_pos'] = 'above'
bz = PulseChannel(label=r'B$_\mathsf{Z}$', color=thc.bc, label_fs=16, **standard_options)

chans = [0]*3
bx_chan = 2;		bz_chan = 1;		gz_chan = 0;
chans[bx_chan] = bx;
chans[bz_chan] = bz;
chans[gz_chan] = gz;

chans4 = [0]*4
bx_chan = 2;		by_chan = 3;		 bz_chan = 1;		gz_chan = 0;
chans4[bx_chan] = bx;
chans4[by_chan] = by;
chans4[bz_chan] = bz;
chans4[gz_chan] = gz;

# Diffusion pulse sequences
dif2 = PulseSequence(chans)						# 2 pulses
dif4 = PulseSequence(chans4)					# 4 pulses

# Have to build the sequence in order. Start with a margin
lmarg_dur = 20;	rmarg_dur = 30
pi_dur = 30
tau_dur = 100

# The labels
rcParams['text.usetex']=False
amsmath = r'\usepackage{xfrac}'
if amsmath not in rcParams['text.latex.preamble']:
	rcParams['text.latex.preamble'].append(amsmath)
label90 = r'$\sfrac{\pi}{2}$' if rcParams['text.usetex'] else r'$\frac{\pi}{2}$'
label180 = r'$\pi$'
labeltauhalf = r'$\sfrac{\tau}{2}$' if rcParams['text.usetex'] else r'$\frac{\tau}{2}$'
labeltau = r'$\tau$'

### 2-pulse version
# The pulses
off_stage = dif2.make_stage(None)
x_stage = dif2.make_stage([bx_chan])
z_stage = dif2.make_stage([bz_chan, gz_chan])
dif2.add_stage(duration=lmarg_dur, states=off_stage)
dif2.add_stage(duration=pi_dur/2, states=x_stage, 
               label=label90, labelchan=bx_chan)
dif2.add_stage(duration=tau_dur/2, states=z_stage,
               label=labeltauhalf, labelchan=bz_chan)
sec_start = dif2.add_stage(duration=pi_dur, states=x_stage,
               label=label180, labelchan=bx_chan)
dif2.add_stage(duration=tau_dur, states=z_stage,
               label=labeltau, labelchan=bz_chan)
dif2.add_stage(duration=pi_dur, states=x_stage,
               label=label180, labelchan=bx_chan)
dif2.add_stage(duration=tau_dur, states=z_stage,
               label=labeltau, labelchan=bz_chan)
sec_end = dif2.add_stage(duration=pi_dur, states=x_stage,
               label=label180, labelchan=bx_chan)
dif2.add_stage(duration=tau_dur/2, states=z_stage,
               label=labeltauhalf, labelchan=bz_chan)
dif2.add_stage(duration=pi_dur/2, states=x_stage, 
               label=label90, labelchan=bx_chan)
dif2.add_stage(duration=rmarg_dur, states=off_stage)

# Repeated sections
rep_sec = dif2.add_section(sec_start, label='N', color=thc.rc, lw=1.5)
dif2.end_section(sec_end, section=rep_sec)

### 4-pulse version
off_stage = dif4.make_stage(None)
x_stage = dif4.make_stage([bx_chan])
y_stage = dif4.make_stage([by_chan])
z_stage = dif4.make_stage([bz_chan, gz_chan])

dif4.add_stage(duration=lmarg_dur, states=off_stage)
dif4.add_stage(duration=pi_dur/2, states=x_stage, 
               label=label90, labelchan=bx_chan)
dif4.add_stage(duration=tau_dur/2, states=z_stage,
               label=labeltauhalf, labelchan=bz_chan)
sec_start = dif4.add_stage(duration=pi_dur, states=x_stage,
               label=label180, labelchan=bx_chan)
dif4.add_stage(duration=tau_dur, states=z_stage,
               label=labeltau, labelchan=bz_chan)
dif4.add_stage(duration=pi_dur, states=x_stage,
               label=label180, labelchan=bx_chan)
dif4.add_stage(duration=tau_dur, states=z_stage,
               label=labeltau, labelchan=bz_chan)
sec_end = dif4.add_stage(duration=pi_dur, states=x_stage,
               label=label180, labelchan=bx_chan)
dif4.add_stage(duration=tau_dur/2, states=z_stage,
               label=labeltauhalf, labelchan=bz_chan)
dif4.add_stage(duration=pi_dur/2, states=x_stage, 
               label=label90, labelchan=bx_chan)
dif4.add_stage(duration=rmarg_dur, states=off_stage)

figsize=(6.5, 3)
dpi = 150
[fig, fnum] = tfs.make_new_fig(fnum, close_old=True, return_fig=True, 
                               figsize=figsize, dpi=dpi)

make_pulse_sequence(fig, dif2,
                    chan_label_weight=0.03,
                    chan_label_right=False)

if SaveFigs:
     tfs.save_figs('../relaxometry/diffusion_sequence_py', formats=['svg', 'eps', 'pdf'])

show()