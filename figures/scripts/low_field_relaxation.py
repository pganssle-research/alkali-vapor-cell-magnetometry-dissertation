"""
Generate figures for low-field relaxation dynamics.

@author Paul J. Ganssle
"""

from scipy.io import loadmat
from scipy.interpolate import interp1d

import theme.colors as thc
import theme.figure_settings as tfs
import os

from matplotlib.pyplot import *
from numpy import squeeze, argsort, append

#
# General setup
#
fnum = 0
fsize=(8, 4)

tfs.setup_environment(figsize=fsize,
                      dpi=150,
                      legend_fs=9,
                      tick_fs=8,
                      label_fs=10,
                      ann_fs=12,
                      n_legend_points=1)

loc = '../relaxometry/'

dcol = thc.rc
wcol = thc.bc

#
# Helper functions and classes
#
class RelaxFieldDataSet:
    def __init__(self, B, Tn, Tne):
        self.B = B*1e3;
        self.Tn = Tn;
        self.Tne = Tne;

        # Derived quantities
        self.f = 4257*self.B
        self.w = self.f*2*pi

        self.R = 1/Tn
        self.Re = self.Tne/(self.Tn**2)       

get_b = lambda otn: squeeze(otn['prog']['aovals'][0,0])*double(otn['disp']['z_cal'])
get_tn = lambda otn, tn, tne: (squeeze(otn['fit'][tn][0,0]),
                               squeeze(otn['fit'][tne][0,0]))
get_t2 = lambda otn: get_tn(otn, 't2', 't2e')
get_t1 = lambda otn: get_tn(otn, 't1', 't1e')

#
# Data loading
#
"""         Water           """
o = loadmat('data_sets/WaterFieldLite.mat')

# Merge the "cold" data sets
ot2f = o['ot2f'][0,0]            # Full frequency
ot2fs = o['ot2fs'][0,0]          # Very low frequency
ot2fl = o['ot2fl'][0,0]          # Longer span

B = get_b(ot2f)
[T2, T2e] = get_t2(ot2f)

Bs = get_b(ot2fs)
[T2s, T2es] = get_t2(ot2fs)

Bl = get_b(ot2fl)
[T2l, T2el] = get_t2(ot2fl)

B = append(B, Bs)
B = append(B, Bl)
T2 = append(T2, T2s)
T2 = append(T2, T2l)
T2e = append(T2e, T2es)
T2e = append(T2e, T2el)

ii = argsort(B)

B = B[ii];      T2 = T2[ii];        T2e = T2e[ii]

T2wc = RelaxFieldDataSet(B, T2, T2e)

# Now the "cold" T1
ot1f = o['ot1f'][0,0]
B = get_b(ot1f)
[T1, T1e] = get_t1(ot1f)

T1wc = RelaxFieldDataSet(B, T1, T1e)

# Now the "hot" T2
ot2fh = o['ot2fh'][0,0]
B = get_b(ot2fh)
[T2, T2e] = get_t2(ot2fh)

T2wh = RelaxFieldDataSet(B, T2, T2e)

"""         Methanol        """
o = loadmat('data_sets/MethanolFieldLite.mat')
# Merge the data sets, as above.
ot2f = o['ot2f'][0,0]            # Full frequency
ot2fs = o['ot2fs'][0,0]          # Very low frequency

B = get_b(ot2f)
[T2, T2e] = get_t2(ot2f)

Bs = get_b(ot2fs)
[T2s, T2es] = get_t2(ot2fs)

B = append(B, Bs)
T2 = append(T2, T2s)
T2e = append(T2e, T2es)

ii = argsort(B)

B = B[ii];      T2 = T2[ii];        T2e = T2e[ii]

T2me = RelaxFieldDataSet(B, T2, T2e)

"""         Ethanol         """
o = loadmat('data_sets/EthanolFieldLite.mat')
ot2f = o['ot2f'][0,0]            
B = get_b(ot2f)
[T2, T2e] = get_t2(ot2f)

T2et = RelaxFieldDataSet(B, T2, T2e)

"""         Octane          """
o = loadmat('data_sets/OctaneFieldLite.mat')
ot2f = o['ot2f'][0,0]           
B = get_b(ot2f)
[T2, T2e] = get_t2(ot2f)

T2oc = RelaxFieldDataSet(B, T2, T2e)

"""         Limonene        """
o = loadmat('data_sets/LimoneneFieldLite.mat')
ot2f = o['ot2f'][0,0]            
B = get_b(ot2f)
[T2, T2e] = get_t2(ot2f)

T2li = RelaxFieldDataSet(B, T2, T2e)

#
#   Plotting
#

err_alpha = 0.75
line_options = {'alpha':err_alpha,
                'ls':'--'}
# Compare "cold" water T1 and T2
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True, figsize=(13, 4))

p2, = plot(T2wc.B, T2wc.R, 'o', ms=9, color=thc.bc)
plot(T2wc.B, T2wc.R, color=thc.bc, lw=1.5, **line_options)

p1, = plot(T1wc.B, T1wc.R, '^', ms=9, color=thc.gc)
plot(T1wc.B, T1wc.R, color=thc.gc, lw=1.5, **line_options)

bmax = max([max(T1wc.B), max(T2wc.B)])
bmax = min([bmax, 600])
xlim(0, bmax)

xlabel('Bias Field (mG)')
ylabel('Relaxation Rate ($s^{-1}$)')

legend([p1, p2], 
       ['Spin-Lattice Relaxation Rate (1/$T_1$)', 'Spin-Spin Relaxation Rate (1/$T_2$)'],
       loc='best')

tight_layout()
tfs.save_figs(os.path.join(loc,'water-relaxation-t1t2'), formats=['pdf', 'eps', 'png'])

# Compare "cold" and "hot" T2
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True, figsize=(6.5, 4))

plot(T2wc.B, T2wc.R, color=thc.bc, **line_options)
pc, = plot(T2wc.B, T2wc.R, 'o', color=thc.bc)

plot(T2wh.B, T2wh.R, color=thc.rc, **line_options)
ph, = plot(T2wh.B, T2wh.R, 's', color=thc.rc)

bmax = max([max(T2wh.B), max(T2wc.B)])
bmax = min([bmax, 750])
xlim(0, bmax)

xlabel('Bias Field (mG)')
ylabel('Spin-Spin Relaxation Rate 1/$T_2$ ($s^{-1}$)')

legend([pc, ph], 
       ['36$^\circ$C', '48$^\circ$C'],
       loc='best')

tight_layout()
tfs.save_figs(os.path.join(loc,'water-relaxation-t2temp'), formats=['pdf', 'eps', 'png'])

# Compare water to various other solvents
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True, figsize=(6.5, 4))

phs = []
legends = []

wcol = thc.bc
mecol = thc.gc
etcol = thc.pc
occol = thc.rc
licol = thc.oc

plot(T2wc.B, T2wc.R, color=wcol, **line_options)
pw, = plot(T2wc.B, T2wc.R, 'o', color=wcol)
legends.append('Water');        phs.append(pw)

plot(T2me.B, T2me.R, color=mecol, **line_options)
pme, = plot(T2me.B, T2me.R, 'v', color=mecol)
legends.append('Methanol');        phs.append(pme)

plot(T2et.B, T2et.R, color=etcol, **line_options)
pet, = plot(T2et.B, T2et.R, '^', color=etcol)
legends.append('Ethanol');        phs.append(pet)

plot(T2oc.B, T2oc.R, color=occol, **line_options)
poc, = plot(T2oc.B, T2oc.R, '8', color=occol)
legends.append('Octane');        phs.append(poc)

plot(T2li.B, T2li.R, color=licol, **line_options)
pli, = plot(T2li.B, T2li.R, 'D', color=licol)
legends.append('Limonene');        phs.append(pli)

xlabel('Bias Field (mG)')
ylabel('Spin-Spin Relaxation Rate 1/$T_2$ ($s^{-1}$)')

bmax = max([max(T2wc.B), max(T2et.B), max(T2me.B), max(T2oc.B), max(T2li.B)])
bmax = min([bmax, 600])
xlim(0, bmax)

legend(phs, legends, loc='best')

tight_layout()
tfs.save_figs(os.path.join(loc,'solvent-relaxation-t2f'), formats=['pdf', 'eps', 'png'])

show()