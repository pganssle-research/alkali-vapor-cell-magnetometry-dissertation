"""
Generate pulse sequence diagrams for demonstration purposes.

@author Paul J. Ganssle
"""
import scipy.io;
from numpy import *;

from matplotlib import rcParams, rc
from matplotlib.pyplot import *;

from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

import matplotlib.gridspec as gridspec

import theme.colors as thc
import theme.figure_settings as tfs

class PulseChannel:
    """
    Group together a few properties of each channel.

    Passing None to any of the attributes uses the global values.
    """

    _label_above = 0
    _label_below = 1
    _label_center = 2
    _label_bottom = 3
    _label_off = 4

    def __init__(self, weight=1,
                 fill_under=False,
                 fill_alpha=1.0,
                 labelweight=0.25,
                 label=None,
                 color=None,
                 lw=None,
                 label_color=None,
                 label_fs='medium',
                 label_fw='regular',
                 label_ha='center',
                 label_va='center',
                 label_pos='above',
                 clabel_va='center',
                 clabel_ha='center',
                 clabel_color=None,
                 clabel_fs='x-large',
                 clabel_fw='regular'):

        label_pos = label_pos.lower()
        __label_pos = {'above':self._label_above,
                       'below':self._label_below,
                       'center':self._label_center,
                       'bottom':self._label_bottom,
                       'off':self._label_off}

        if label_pos not in __label_pos.keys():
            raise ValueError('Invalid label position. '+\
                             'Label position must be one of these keys: '+\
                             ','.join(__label_pos.keys()))
        else:
            label_pos = __label_pos[label_pos]

        if label_color is None and color is not None:
            label_color = color
        
        if clabel_color is None and color is not None:
            clabel_color = color

        self.label = label
        self.color = color     
        self.label_color = label_color
        self.clabel_color = clabel_color
        
        self.lw=lw
        self.weight=weight
        self.labelweight = labelweight
        self.fill_under = fill_under
        self.fill_alpha = fill_alpha

        self.label_fs = label_fs
        self.label_fw = label_fw
        self.label_ha = label_ha
        self.label_va = label_va
        self.label_pos = label_pos

        self.clabel_fs = clabel_fs
        self.clabel_fw = clabel_fw
        self.clabel_ha = clabel_ha
        self.clabel_va = clabel_va

        self.ax = None

    def has_label(self):
        """
        Check if the channel has a label or not.
        """
        return self.label is not None

class PulseSequence:
    """
    Pulse sequence container class.
    """
    class Pulse:
        """
        Simple class to contain information about the positioning of each pulse
        """
        def __init__(self, start, width=None, end=None):
            if width is None and end is None:
                raise ValueError('Must provide either a width or the end location')

            self.start = start
            if width is not None:
                self.width = width
                self.end = start+width
            else:
                self.end = end
                self.width = end-start

            self.center = mean([self.start, self.end])

    class Section:
        """
        Simple class for containing repeats
        """
        def __init__(self, start,
                     end=None,
                     label=None,
                     color=thc.rc,
                     lw=1.0,
                     chan_bottom=None,
                     chan_top=None,
                     offset_before=0.1,
                     offset_after=0.1,
                     offset_above=0.03,
                     offset_below=0.03,
                     text_offset_right=0.02,
                     text_offset_below=0.01,
                     lip_frac=0.03,):


            self.start = start
            self.end  = end
            self.label = label
            self.offset_before = offset_before
            self.offset_after = offset_after
            self.offset_above = offset_above
            self.offset_below = offset_below
            self.text_offset_right = text_offset_right
            self.text_offset_below = text_offset_below
            self.lip_frac = lip_frac
            self.chan_bottom=chan_bottom
            self.chan_top=chan_top
            self.color=color
            self.lw=lw

        def add_end(self, end):
            self.end=end

    class Line:
        """
        Simple class for single events
        """
        def __init__(self, pos,
                     label=None,
                     color=thc.rc,
                     lw=1.0,
                     label_fs=10,
                     label_fw='regular',
                     label_ha = 'center',
                     label_va = 'top',
                     chan_bottom=None,
                     chan_top=None,
                     offset_before=0.00,
                     offset_above=0.03,
                     offset_below=0.03,
                     text_offset_right=0.0,
                     text_offset_below=0.01):

            self.pos = pos
            self.label = label
            self.label_fs = label_fs
            self.label_fw = label_fw
            self.label_ha = label_ha
            self.label_va = label_va

            self.offset_before=offset_before
            self.offset_above = offset_above
            self.offset_below = offset_below
            
            self.text_offset_right = text_offset_right
            self.text_offset_below = text_offset_below
            
            self.chan_bottom=chan_bottom
            self.chan_top=chan_top
            
            self.color=color
            self.lw=lw

    def __init__(self, channels):
        self.channels = channels
        self._num_chans = len(channels)
        self.chan_labels = [False]*self._num_chans
        self.up_to_date=False

        # Pulse details
        self._num_stages = 0
        self.durations = []                                     # Pulse durations
        self.stages = [[] for ii in range(0, self._num_chans)]  # List of on/off states.

        # Plotting details
        self._total_duration = 0.0
        self.pulse_widths = []
        self.pulses=[]                      # Each entry should be (start, end, center)

        self.labels = []                                        # Label should be None or a string
        self.labelchans = []                                    # Which channel is labeled
        self.label_fs =[]

        # Repeated sections
        self.sections = []
        self.sections_start = []            # Keep track of what sections start here.
        self.sections_end = []

        # Lines
        self.lines = []

    def add_stage(self, duration, states,
        label=None, labelchan=None, label_fs=None,
        update=True):
        """
        Add a segment to the pulse sequence.
        """
        # Input validation
        if duration <= 0:
            raise ValueError('Duration must be greater than zero. You have supplied: '+\
                             str(duration))

        if len(states) < self._num_chans:
            raise ValueError(('Valid state must be supplied for {:0.0f} states. You have '+
                             'supplied only {:0.0f} states.').format(self._num_chans, len(states)))

        if label is not None:
            label = str(label)

        if labelchan is not None:
            self.chan_labels[labelchan] = True

        self._num_stages += 1
        self.labels.append(label)
        self.label_fs.append(label_fs)
        self.labelchans.append(labelchan)
        self.durations.append(double(duration))
        self.sections_start.append([])
        self.sections_end.append([])
        for ii, state in enumerate(self.stages):
            state.append(states[ii])

        self.up_to_date = False
        if update:
            self._update_widths()

        return self._num_stages

    def add_line(self,
                    stage_num,
                    label=None,
                    color=None,
                    lw=None,
                    label_fs=None,
                    label_fw=None,
                    label_ha=None,
                    label_va=None,
                    offset_before=None,
                    offset_above=None,
                    offset_below=None,
                    text_offset_right=None,
                    text_offset_below=None,
                    chan_bottom=None,
                    chan_top=None):
        """
        Adds a new repeated subsequence.
        """
        # Create a kwargs dict to pass on args, but don't pass on None values,
        # as those might overwrite the defaults coded into the constructor prototype
        options_dict = {
            'label':label,
            'color':color,
            'lw':lw,
            'label_fs':label_fs,
            'label_fw':label_fw,
            'label_ha':label_ha,
            'label_va':label_va,
            'offset_before':offset_before,
            'offset_below':offset_below,
            'offset_above':offset_above,
            'text_offset_right':text_offset_right,
            'text_offset_below':text_offset_below,
            'chan_bottom':chan_bottom,
            'chan_top':chan_top,
        }

        for key in options_dict.keys():
            if options_dict[key] is None:
                del options_dict[key]

        cline = self.Line(stage_num, **options_dict)
        self.lines.append(cline)
        
        return cline

    def add_section(self,
                    stage_num,
                    label=None,
                    color=None,
                    lw=None,
                    offset_before=None,
                    offset_after=None,
                    offset_above=None,
                    offset_below=None,
                    text_offset_right=None,
                    text_offset_below=None,
                    lip_frac=None,
                    chan_bottom=None,
                    chan_top=None):
        """
        Adds a new repeated subsequence.
        """
        # Create a kwargs dict to pass on args, but don't pass on None values,
        # as those might overwrite the defaults coded into the constructor prototype
        options_dict = {
            'label':label,
            'color':color,
            'offset_before':offset_before,
            'offset_after':offset_after,
            'offset_below':offset_below,
            'offset_above':offset_above,
            'text_offset_right':text_offset_right,
            'text_offset_below':text_offset_below,
            'lip_frac':lip_frac,
            'chan_bottom':chan_bottom,
            'chan_top':chan_top,
            'lw':lw
        }

        for key in options_dict.keys():
            if options_dict[key] is None:
                del options_dict[key]

        section = self.Section(stage_num, **options_dict)
        self.sections.append(section)
        self.sections_start[stage_num].append(section)

        return section


    def end_section(self, stage_num, section=None):
        """
        Adds the end of the section to one of the sections in the sequence.

        If section is a number, then it ends the Nth section. Otherwise it
        finds the last unclosed section and closes it. If a Section is passed, then it ends that
        one.
        """

        if len(self.sections) < 1:
            raise Exception('There are no sections to end.')

        if section is None:
            # Find the most recent unopened section. It makes sense to iterate like this:
            #           { S1, S2, S3, [Z1, Z2, Z3]*N, S4, S5}*M
            # But it makes no sense to iterate like this:
            #           { S1, S2, S3, [Z1, Z2, Z3}*M, S4, S5]*N
            # So we'll assume it's FIFO
            section_found = False
            for ii in reversed(range(0, len(self.sections))):
                if self.sections[ii].end is None:
                    # Found an unclosed one
                    section_found = True
                    section = self.sections[ii]
                    break

            if not section_found:
                raise Exception('All sections are already closed.')

        section.end = stage_num
        self.sections_end[stage_num].append(section)

    def _update_widths(self):
        """
        Recalculate the width of each pulse, as a fraction of the total length of the sequence
        """

        self._total_duration = double(sum(self.durations))
        self.pulse_widths = []
        self.pulses = []

        cend = 0.0
        for dur in self.durations:
            cstart = cend                       # The start of this pulse is the end of the last one
            width = dur/self._total_duration    # Width of the pulse is a fraction of total duration

            pulse = self.Pulse(cstart, width=width)
            cend = pulse.end

            self.pulse_widths.append(width)
            self.pulses.append(pulse)

        self._up_to_date = True

    def make_stage(self, chans_on):
        """
        Creates a bool array with an entry for each channel, false if the channel index (0-based)
        is not in the list chans_on
        """
        if chans_on is None:
            chans_on = []

        if isinstance(chans_on, int):
            chans_on = [chans_on]

        return [i in chans_on for i in range(0, self._num_chans)]

def make_pulse_sequence(fig, seq,
                        show_labels=True,
                        show_chan_labels=True,
                        chan_label_weight=0.2,
                        chan_label_right=True,
                        tight_pad=0.5,
                        hspace=0.5,
                        wspace=0.05):
    """
    Makes the pulse sequence, given a sequence seq
    """

    # Store the axes.
    caxs = [None]*seq._num_chans;       laxs = [None]*seq._num_chans
    nlabels = sum([1 if x else 0 for x in seq.chan_labels])

    nlabels = 0
    column_on = [False for x in range(0, 2*seq._num_chans+1)]
    label_for = [[] for x in range(0, 2*seq._num_chans+1)]
    label_row = [None for x in range(0, seq._num_chans)]
    is_channel = [None for x in range(0, 2*seq._num_chans+1)]

    for ii, chan in enumerate(seq.channels):
        if show_labels and seq.chan_labels[ii]:
            if chan.label_pos == chan._label_bottom:
                column_on[-1] = True
                label_row[ii] = len(column_on)-1
                label_for[-1].append(ii)
            elif chan.label_pos == chan._label_above:
                column_on[2*ii] = True
                label_for[2*ii].append(ii)
                label_row[ii] = 2*ii
            elif chan.label_pos == chan._label_below:
                column_on[2*ii+2] = True
                label_for[2*ii+2].append(ii)
                label_row[ii] = 2*ii+2
            elif chan.label_pos == chan._label_center:
                label_row[ii] = 2*ii+1
                label_for[2*ii+1].append(ii)

        column_on[2*ii+1] = True        # Always true.
        is_channel[2*ii+1] = ii

    # Delete the columns we're not using
    columns_on_buff = []
    label_for_buff = []
    is_channel_buff = []
    for ii in range(0, len(column_on)):
        if not column_on[ii]:
            for jj in range(ii, len(label_for)):
                for lc in label_for[jj]:
                    label_row[lc] -= 1
        else:
            columns_on_buff.append(True)
            label_for_buff.append(label_for[ii])
            is_channel_buff.append(is_channel[ii])

    column_on = columns_on_buff
    label_for = label_for_buff
    is_channel = is_channel_buff

    lrows = label_row
    crows = [None for x in range(0, seq._num_chans)]
    for ii in range(0, len(is_channel)):
        if is_channel[ii] is not None:
            crows[is_channel[ii]] = ii


    # Build the gridspec
    ngsrows = sum([1 if x else 0 for x in column_on])
    ngscols = 2 if show_chan_labels else 1
    row_weights = [0.0 for x in column_on]

    for ii in range(0, len(label_for)):
        if column_on[ii]:
            if is_channel[ii] is not None:
                row_weights[ii] = seq.channels[is_channel[ii]].weight
            elif len(label_for[ii]) > 0:
                label_weights = []
                for c in label_for[ii]:
                    label_weights.append(seq.channels[c].labelweight)

                row_weights[ii] = max(label_weights)
            
    if show_chan_labels:
        ngscols = 2
        col_weights = [0, 0]
        (chan_col, clab_col) = (0, 1) if chan_label_right else (1, 0)

        col_weights[chan_col] = 1.0
        col_weights[clab_col] = chan_label_weight
    else:
        ngscols = 1
        col_weights = [1.0]
        chan_col = 0

    gs = gridspec.GridSpec(ngsrows, ngscols, width_ratios=col_weights, height_ratios=row_weights)

    pulse_starts = [pulse.start for pulse in seq.pulses]
    pulse_ends = [pulse.end for pulse in seq.pulses]


    axs = [None for x in range(0, len(row_weights))]

    for ii in range(0, seq._num_chans):
        crow = crows[ii]
        chan = seq.channels[ii]

        # Make the channel
        if axs[crow] is None:
            cax = subplot(gs[crow, chan_col])
            axs[crow] = cax
        else:
            cax = axs[crow]

        caxs[ii] = cax
        
        cax.axis('off')

        # Draw a single line that follows the pulse correctly
        cstate = False
        xdata = [0];            ydata = [0]             # Starts at 0,0
        xfdata = [0];           yfdata = [0]

        for jj in range(0, seq._num_stages):
            cstage = seq.stages[ii][jj]
            if cstage != cstate:
                # This is a transition, so add a new line segment
                xdata.append(seq.pulses[jj].start)
                ydata.append(1.0 if cstage else 0.0)

                xfdata.append(seq.pulses[jj].start)
                yfdata.append(1.0 if cstate else 0.0)

                xfdata.append(seq.pulses[jj].start)
                yfdata.append(1.0 if cstage else 0.0)
                cstate = cstage

        xfdata.append(1)
        yfdata.append(1.0 if cstage else 0.0)

        xdata.append(1)             # The very end is always 1.0
        ydata.append(1.0 if cstage else 0.0)

        if chan.fill_under:
            fill_between(xfdata, yfdata, color=chan.color, alpha=chan.fill_alpha)
            plot([min(xfdata), max(xfdata)], [0,0], lw=chan.lw, color=chan.color, 
                 transform=cax.transAxes, clip_on=False, alpha=chan.fill_alpha)


        line = Line2D(xdata, ydata, c=chan.color,
                      drawstyle='steps-post',
                      solid_joinstyle='round',
                      lw=chan.lw,
                      clip_on=False,
                      transform=cax.transAxes)
        cax.add_line(line)
        
        xlim([0, 1])
        ylim([0, 1])

        # Add in the labels now, if applicable
        if show_labels and seq.chan_labels[ii]:
            lrow = lrows[ii]

            if axs[lrow] is not None:
                lax = axs[lrow]
            else:
                lax = subplot(gs[lrow, chan_col])
                lax.axis('off')
        
            laxs[ii] = lax

            for jj in range(0, seq._num_stages):
                label = seq.labels[jj]
                if label is None or seq.labelchans[jj] != ii:
                    continue
                
                if seq.label_fs[jj] is not None:
                    label_fs = seq.label_fs[jj]
                else:
                    label_fs = chan.label_fs

                lax.text(seq.pulses[jj].center, 0.5,
                        label,
                        va=chan.label_va,
                        ha=chan.label_ha,
                        clip_on=False,
                        color=chan.label_color,
                        fontsize=label_fs,
                        fontweight=chan.label_fw,
                        transform=lax.transAxes)

        if show_chan_labels and chan.has_label():
            clax = subplot(gs[crow, clab_col])
            clax.axis('off')

            clax.text(0.5, 0.5,
                      chan.label,
                      va=chan.clabel_va,
                      ha=chan.clabel_ha,
                      clip_on=False,
                      color=chan.clabel_color,
                      fontsize=chan.clabel_fs,
                      fontweight=chan.clabel_fw,
                      transform=clax.transAxes)

    # Need to do this before working in transFigure axes
    tight_layout(tight_pad)
    subplots_adjust(hspace=hspace, wspace=wspace)
    
    '''
    # If we're showing the labels below the chart, we'll do it here
    if show_labels_below:
        lax = subplot(gs[-1, chan_col])
        lax.axis('off')

        for ii in range(0, seq._num_stages):
            label = seq.labels[ii]
            if label is None:
                continue

            lax.text(seq.pulses[ii].center, 0.5,
                     label,
                     va='center',
                     ha='center',
                     clip_on=False,
                     color=seq.channels[seq.labelchans[ii]].color,
                     fontsize='medium',
                     fontweight='regular',
                     transform=lax.transAxes)
    '''

    # Plot the repeated sections portion of this.
    c_width_tf = caxs[-1].get_position().width          # Width of channels in transFigure units
    c_left_tf = caxs[-1].get_position().xmin            # Left side of channel in tF units
    c_bot_tf = caxs[-1].get_position().ymin
    c_top_tf = caxs[0].get_position().ymax              # Top of the highest channel, in tF units

    for ii, section in enumerate(seq.sections):
        # Do this in transFigure units, to span all the channels.
        start_stage = section.start-1
        end_stage = section.end-1

        ss = seq.pulses[start_stage]
        es = seq.pulses[end_stage]

        l_ta = ss.end-ss.width*section.offset_before
        r_ta = es.start+es.width*section.offset_after

        l_tf = c_left_tf + l_ta*c_width_tf          # Transform to transFigure axes
        r_tf = c_left_tf + r_ta*c_width_tf

        # Draw the lines now
        lip_frac = section.lip_frac
        offset_below = section.offset_below
        offset_above = section.offset_above

        text_offset_right = section.text_offset_right
        text_offset_below = section.text_offset_below

        sec_wid = r_tf-l_tf
        sec_ho = c_top_tf-c_bot_tf

        lip_w = (r_tf-l_tf)*lip_frac
        t_tf = c_top_tf+sec_ho*offset_above
        b_tf = c_bot_tf-sec_ho*offset_below
        sec_h = t_tf-b_tf

        leftLine = Line2D([l_tf+lip_w, l_tf, l_tf, l_tf+lip_w],
                          [b_tf, b_tf, t_tf, t_tf],
                          solid_joinstyle='round',
                          clip_on=False,
                          color=section.color,
                          lw=section.lw,
                          transform=fig.transFigure)

        rightLine = Line2D([r_tf-lip_w, r_tf, r_tf, r_tf-lip_w],
                          [b_tf, b_tf, t_tf, t_tf],
                          solid_joinstyle='round',
                          clip_on=False,
                          color=section.color,
                          lw=section.lw,
                          transform=fig.transFigure)

        # Add these to an arbitrary axis, I guess.
        caxs[-1].add_line(leftLine)
        caxs[-1].add_line(rightLine)

        # Add the label
        if section.label is not None:
            caxs[-1].text(r_tf+text_offset_right*sec_wid,
                 b_tf-text_offset_below*sec_ho,
                 section.label,
                 va='center',
                 ha='left',
                 fontsize='small',
                 color=section.color,
                 clip_on=False,
                 transform=fig.transFigure)

    # Add any single-event lines
    c_width_tf = caxs[-1].get_position().width          # Width of channels in transFigure units
    c_left_tf = caxs[-1].get_position().xmin            # Left side of channel in tF units
    c_bot_tf = caxs[-1].get_position().ymin
    c_top_tf = caxs[0].get_position().ymax              # Top of the highest channel, in tF units

    for ii, cline in enumerate(seq.lines):
        # Do this in transFigure units, to span all the channels
        if cline.pos >= len(seq.pulses):
            pos = seq.pulses[-1].end
        else:
            pos = seq.pulses[cline.pos].start

        # Draw the lines now
        offset_below = cline.offset_below
        offset_above = cline.offset_above
        offset_before = cline.offset_before
        
        text_offset_right = cline.text_offset_right
        text_offset_below = cline.text_offset_below

        
        theLine = Line2D([c_left_tf+(pos-offset_before)*c_width_tf]*2,
                          [c_bot_tf-offset_below, c_top_tf+offset_above],
                          solid_capstyle='projecting',
                          clip_on=False,
                          color=cline.color,
                          lw=cline.lw,
                          transform=fig.transFigure)

        # Add to an arbitrary axis
        caxs[-1].add_line(theLine)

        # Add the label
        if cline.label is not None:
            caxs[-1].text(c_left_tf+(pos+text_offset_right)*c_width_tf,
                 c_bot_tf-text_offset_below,
                 cline.label,
                 va=cline.label_va,
                 ha=cline.label_ha,
                 fontsize=cline.label_fs,
                 fontweight=cline.label_fw,
                 color=cline.color,
                 clip_on=False,
                 transform=fig.transFigure)



