"""
Generate figures from the decane/water T1, T2 and diffusion data

@author Paul J. Ganssle
"""

from scipy.io import loadmat

import theme.colors as thc
import theme.figure_settings as tfs
import os

from matplotlib.pyplot import *
from numpy import squeeze

# Load data
o = loadmat('data_sets/DecaneWaterRelaxation2.mat')
lsfdt1 = o['lsfdt1'];		lsfdt2 = o['lsfdt2'];
lsfwt1 = o['lsfwt1'];		lsfwt2 = o['lsfwt2'];
lsfdwt1 = o['lsfdwt1'];		lsfdwt2 = o['lsfdwt2'];
lsfdd = o['lsfdd'];         lsfwd = o['lsfwd']
lsfdwd = o['lsfdwd']

# Laplace inversions
o1dt1 = o['o1dt1'];			o1dt2 = o['o1dt2']
o1wt1 = o['o1wt1'];			o1wt2 = o['o1wt2']
o1dwt1=o['o1dwt1'];			o1dwt2 = o['o1dwt2']

# T1 - Fit Parameters
t1d = lsfdt1['tau'][0][0][0][0];            t1dA = lsfdt1['tau'][0][0][0][1];
t1w = lsfwt1['tau'][0][0][0][0];            t1wA = lsfwt1['tau'][0][0][0][1];
t1dw = lsfdwt1['tau'][0][0][0][0::2];       t1dwA = lsfdwt1['tau'][0][0][0][1::2];

# T1 - Fit Data and curves
td1 = squeeze(lsfdt1['t'][0][0]);           
tw1 = squeeze(lsfwt1['t'][0][0]);           
tdw1 = squeeze(lsfdwt1['t'][0][0]);         

cd1 = squeeze(lsfdt1['c'][0][0]);           cfd1 = squeeze(lsfdt1['cf'][0][0]);
cw1 = squeeze(lsfwt1['c'][0][0]);           cfw1 = squeeze(lsfwt1['cf'][0][0]);
cdw1 = squeeze(lsfdwt1['c'][0][0]);         cfdw1 = squeeze(lsfdwt1['cf'][0][0]);

#T2 - Fit Parameters
t2d = lsfdt2['tau'][0][0][0][0];            t1dA = lsfdt2['tau'][0][0][0][1];
t2w = lsfwt2['tau'][0][0][0][0];            t1wA = lsfwt2['tau'][0][0][0][1];
t2dw = lsfdwt2['tau'][0][0][0][0::2];       t1dwA = lsfdwt2['tau'][0][0][0][1::2];

# T1 Fit Data and Curves
td2 = squeeze(lsfdt2['t'][0][0]);
tw2 = squeeze(lsfwt2['t'][0][0]);
tdw2 = squeeze(lsfdwt2['t'][0][0]);

cd2 = squeeze(lsfdt2['c'][0][0]);           cfd2 = squeeze(lsfdt2['cf'][0][0]);
cw2 = squeeze(lsfwt2['c'][0][0]);           cfw2 = squeeze(lsfwt2['cf'][0][0]);
cdw2 = squeeze(lsfdwt2['c'][0][0]);         cfdw2 = squeeze(lsfdwt2['cf'][0][0]);

# Diffusion - Fit Parameters
Dd = lsfdd['tau'][0][0][0][0]*1e5;          DdA = lsfdd['tau'][0][0][0][1];
Dw = lsfwd['tau'][0][0][0][0]*1e5;          DwA = lsfwd['tau'][0][0][0][1];
Ddw = lsfdwd['tau'][0][0][0][0::2]*1e5;     DdwA = lsfdwd['tau'][0][0][0][1::2];

# Diffusion - Fit Data and curves
tdd = squeeze(lsfdd['t'][0][0]);    
twd = squeeze(lsfwd['t'][0][0]);
tdwd = squeeze(lsfdwd['t'][0][0]);  

cdd = squeeze(lsfdd['c'][0][0]);            cfdd = squeeze(lsfdd['cf'][0][0]);
cwd = squeeze(lsfwd['c'][0][0]);            cfwd = squeeze(lsfwd['cf'][0][0]);
cdwd = squeeze(lsfdwd['c'][0][0]);          cfdwd = squeeze(lsfdwd['cf'][0][0]);

'''
#This sort of thing was in the original for whatever reason. 
selv = (tw2 <= 12)
tw2 = tw2[selv]
cw2 = cw2[selv]/max(cw2[selv])
cfw2 = cfw2[selv]/max(cw2[selv])
'''

# Colors
bluec='#3953a4';
redc = '#b5130c';
greenc = '#007935';
purplec = '#880085';

bluec = thc.bc
redc = thc.rc
greenc = thc.gc
purplec = thc.pc

# Figure setup
fnum = 0
fsize=(6.5, 4)
tfs.setup_environment(figsize=fsize,
                      dpi=150,
                      legend_fs=9,
                      tick_fs=7,
                      n_legend_points=1)


lb1c = 'Fit: ${:0.2f}e^{{-t/{:0.2f}}}$'
lb2c = 'Fit: ${:0.2f}e^{{-t/{:0.2f}}}+{:0.2f}e^{{-t/{:0.2f}}}$'

lb1c = 'Fit: $e^{{-t/{:0.2f}}}$'
lb2c = 'Fit: $e^{{-t/{:0.2f}}}+e^{{-t/{:0.2f}}}$'

d_m = 'o';		d_l = '-';
w_m = 's';		w_l = '-';
dw_m = '^';		dw_l = '-';

dcol = redc
wcol = bluec
dwcol = purplec

dzorder = 3
wzorder = 2
dwzorder = 1

ms = 8
lw = 2

loc = '../relaxometry'

# Decane/Water T1
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True)
clf()
[pdf, pd] = plot(td1, cfd1/max(cd1), d_l, 
                 td1, cd1/max(cd1), d_m, 
                 color=dcol, zorder=dzorder, ms=ms, lw=lw);
[pwf, pw] = plot(tw1, cfw1/max(cw1), w_l,
                 tw1, cw1/max(cw1), w_m,  
                 color=wcol, zorder=wzorder, ms=ms, lw=lw);
[pdwf, pdw] = plot(tdw1, cfdw1/max(cdw1), dw_l,
                   tdw1, cdw1/max(cdw1), dw_m, 
                   color=dwcol, zorder=dwzorder, ms=ms, lw=lw);
gca().set_alpha(0);
gcf().set_alpha(0);
xlim(0, 12)
ylim(-0.025, 1.025)
xlabel('Time (s)')
ylabel('Norm. Signal')
#subplots_adjust(left=0.1, bottom=0.16, right=0.97, top=0.96);
tight_layout()

legend([pd, pdf, pw, pwf, pdw, pdwf], 
       ['Decane', lb1c.format(t1d), 
        'Water', lb1c.format(t1w),
        'Decane/Water Mix', lb2c.format(t1dw[0], t1dw[1])],
        loc='best')

'''
legend([pd, pdf, pw, pwf, pdw, pdwf], 
       ['Decane', lb1c.format(t1dA, t1d), 
        'Water', lb1c.format(t1wA, t1w),
        'Decane/Water Mix', lb2c.format(t1dwA[0], t1dw[0], t1dwA[1], t1dw[1])],
        loc='best')
'''

tfs.save_figs(os.path.join(loc,'decane-water-t1'), formats=['pdf', 'eps', 'png'])

# Decane/Water T2
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True)
clf()
[pdf, pd] = plot(td2, cfd2/max(cd2), d_l, 
                 td2, cd2/max(cd2), d_m, 
                 color=dcol, zorder=dzorder, ms=ms, lw=lw);
[pwf, pw] = plot(tw2, cfw2/max(cw2), w_l, 
                 tw2, cw2/max(cw2), w_m, 
                 color=bluec, zorder=wzorder, ms=ms, lw=lw);
[pdwf, pdw] = plot(tdw2, cfdw2/max(cdw2), dw_l, 
                   tdw2, cdw2/max(cdw2), dw_m, 
                   color=purplec, zorder=dwzorder, ms=ms, lw=lw);

gca().set_alpha(0);
gcf().set_alpha(0);
ylim(-0.025, 1.025)
xlim(0, 12)
xlabel('Time (s)')
ylabel('Norm. Signal')
#subplots_adjust(left=0.1, bottom=0.16, right=0.97, top=0.96);
tight_layout()

legend([pd, pdf, pw, pwf, pdw, pdwf], 
       ['Decane', lb1c.format(t2d), 
        'Water', lb1c.format(t2w),
        'Decane/Water Mix', lb2c.format(t1dw[0], t2dw[1])],
        loc='best')

'''
legend([pd, pdf, pw, pwf, pdw, pdwf], 
       ['Decane', lb1c.format(t1dA, t1d), 
        'Water', lb1c.format(t1wA, t1w),
        'Decane/Water Mix', lb2c.format(t1dwA[0], t1dw[0], t1dwA[1], t1dw[1])],
        loc='best')
'''
tfs.save_figs(os.path.join(loc,'decane-water-t2'), formats=['pdf', 'eps', 'png'])

# Decane/Water Diffusion
[fig, fnum] = tfs.make_new_fig(fnum, return_fig=True, figsize=(13, 6))

[pdf, pd] = plot(tdd, cfdd/max(cdd), d_l, 
                 tdd, cdd/max(cdd), d_m, 
                 color=dcol, zorder=dzorder, ms=ms, lw=lw);
[pwf, pw] = plot(twd, cfwd/max(cwd), w_l, 
                 twd, cwd/max(cwd), w_m, 
                 color=wcol, zorder=wzorder, ms=ms, lw=lw);
[pdwf, pdw] = plot(tdwd, cfdwd/max(cdwd), dw_l, 
                   tdwd, cdwd/max(cdwd), dw_m, 
                   color=dwcol, zorder=dwzorder, ms=ms, lw=lw);

gca().set_alpha(0);
gcf().set_alpha(0);
ylim(-0.025, 1.025)
xlim(0, max(twd))

lb1c = r'$e^{{-(\gamma G\tau)^2{:0.2f}\cdot 10^{{-5}}\tau/3}}$'
lb2c = r'${:0.2f}e^{{-(\gamma G\tau)^2{:0.2f}\cdot 10^{{-5}}\tau/3}} '+\
       r'+ {:0.2f}e^{{-(\gamma G\tau)^2{:0.2f}\cdot 10^{{-5}}\tau/3}}$'

lb1c = 'D$_{{{}}}$ = {:0.2f}$\cdot$10$^{{-5}} cm^{{2}}\!/s$'
lb2c = 'D$_{{d}}$ = {:0.2f}$\cdot$10$^{{-5}} cm^{{2}}\!/s$, '+\
       'A$_{{d}}$ = {:0.2f}\nD$_{{w}}$ = {:0.2f}$\cdot$10$^{{-5}} cm^{{2}}\!/s$, '+\
       'A$_{{w}}$ = {:0.2f}'

'''
legend([pd, pdf, pw, pwf, pdw, pdwf], 
       ['Decane', lb1c.format(Dd), 
        'Water', lb1c.format(Dw),
        'Decane/Water Mix', lb2c.format(DdwA[0], Ddw[0], DdwA[1], Ddw[1])],
        loc='best')
'''                                    
legend([pd, pdf, pw, pwf, pdw, pdwf], 
       ['Decane', lb1c.format('d', Dd), 
        'Water', lb1c.format('w', Dw),
        'Decane/Water Mix', lb2c.format(Ddw[0], DdwA[0], Ddw[1], DdwA[1])], 
        loc='best')
                                       
xlabel('Gradient Strength (G/cm)')
ylabel('Norm. Signal')
#subplots_adjust(left=0.1, bottom=0.16, right=0.97, top=0.96);

tight_layout()                         

tfs.save_figs(os.path.join(loc,'decane-water-diff'), formats=['pdf', 'eps', 'png'])
					   
show();

