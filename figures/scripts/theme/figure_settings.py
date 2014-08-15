from matplotlib.pyplot import fignum_exists, clf, figure, close
from matplotlib.pyplot import rc
from matplotlib.pyplot import savefig
from matplotlib.pyplot import get_fignums
import theme.colors as thc
import re
import os

# Methods
def make_new_fig(fnum=None, close_old=True, return_fig=False, 
                 use_max_fnum=False, **kwargs):
    """
    Make a new figure, clearing the old one.
    
    Returns the figure after the one you've created.
    """

    if fnum is None:
        fignums = get_fignums()
        max_fnum = max(fignums) if len(fignums) > 0 else 0

        if use_max_fnum:
            fnum = max_fnum+1
        else:
            for fnum in range(0, max_fnum+2):
                if fnum not in fignums:
                    break

    if fignum_exists(fnum) and close_old:
        close(fnum)            # Close the figure if it already exists.

    fig = figure(fnum, **kwargs)
    clf()

    if return_fig:
        return [fig, fnum+1]
    else:
        return fnum+1

def close_all_figs(max_fnum=100):
    """
    Close all the existing figures.
    """

    fnums = get_fignums()

    for fnum in fnums:
        close(fnum)


def save_figs(loc, formats='pdf', **kwargs):
    """
    Saves the file at 'loc' with all the formats iterated in "formats"
    """
    # Check if the path exists and make it if it doesn't.
    basepath = '/'.join(re.split(r'[\\/]',loc)[0:-1])
    if len(basepath) > 0 and not os.path.exists(basepath):
        os.makedirs(basepath)

    # Can pass a string or a list/tuple of strings
    if isinstance(formats, str):
        formats = [formats]

    for cformat in formats:
        savefig(loc+'.'+cformat, format=cformat, **kwargs)


def setup_environment(font_type='sans', 
                      tex_regular=True,
                      dpi=150,
                      figsize=(6,4),
                      title_fs=12,
                      label_fs=8,
                      ann_fs=10,
                      tick_fs=6,
                      legend_fs=12,
                      n_legend_points=3,
                      legend_shadow=False,
                      color_cyc=thc.extended_color_cycle):
    """
    Sets up the environment with the right rcParams and such

    font_type can be any in ['sans', 'sans-serif', 'serif', 'mono']
    Default is 'sans', sans and sans-serif are aliases
    """

    sans_type = 'sans-serif';        serif_type = 'serif';        mono_type = 'monospace';
    font_types = {'sans':sans_type, 'sans-serif':sans_type, 'serif':serif_type, 'mono':mono_type}
    if font_type not in font_types.keys():
        raise KeyError("Invalid font type: "+repr(font_type)+". Valid types are "+
                       ', '.join(["'"+ft+"'" for ft in font_types.keys()]))

    # Figure settings
    cdpi = dpi;                    cfigsize = figsize;
    cfacecolor = 'none';        cedgecolor = 'none';    
    rc('figure', figsize=cfigsize,
                dpi=cdpi,
                facecolor=cfacecolor,
                edgecolor=cedgecolor)
    rc('savefig', dpi=cdpi,
                facecolor=cfacecolor,
                edgecolor=cedgecolor,
                pad_inches=0.05,
                frameon=False,
                bbox='tight')

    # Plot settings
    rc('lines', dash_joinstyle='miter', 
                dash_capstyle='butt', 
                solid_capstyle='butt')

    rc('axes', color_cycle=color_cyc, 
               facecolor='none',
               edgecolor='k',
               titlesize=title_fs,
               labelsize=label_fs,
               labelweight='normal')

    rc('legend', fontsize=legend_fs,
                 numpoints=n_legend_points,
                 shadow=legend_shadow)

    rc('xtick', labelsize=tick_fs)
    rc('ytick', labelsize=tick_fs)

    # Font settings
    ss_fonts = ['Myriad Pro', 'DejaVu Sans', 'Bitstream Vera Sans', 'Helvetica', 'sans-serif']
    s_fonts = ['DejaVu', 'Bitstream Vera Serif', 'New Century Schoolbook L', 'serif']
    m_fonts = ['DejaVu Mono', 'Bistream Vera Mono', 'Andale Mono', 
                'Numbus Mono L', 'Terminal', 'monospace']

    cfont_type = font_types[font_type]
    font_dict = {'family':cfont_type, sans_type:ss_fonts, serif_type:s_fonts, mono_type:m_fonts}
    
    rc('text', usetex=False)
    rc('font', **font_dict)

    if cfont_type == sans_type:
        rc('mathtext', fontset='stixsans', default='regular' if tex_regular else 'sf')
    elif cfont_type == serif_type:
        rc('mathtext', fontset='stix', default='regular' if tex_regular else 'rm')
    elif cfont_type == mono_type:
        rc('mathtext', fontset='stix', default='regular' if tex_regular else 'tt')
