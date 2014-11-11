#!/usr/bin/env python
# By Thijs Coenen (2012-2013) for PhD thesis research at the Universiteit van
# Amsterdam.
'''
Condensed single pulse plots from PRESTO .singlepulse/.inf files.

Second large iteration, meant to produce plots that can be included in thesis,
so no black backgrounds and glitches fixed.
'''
# Python standard library imports
from __future__ import division
import os
import sys
import optparse
import StringIO
from base64 import encodestring
import re
from collections import defaultdict

# Standard 3rd party libraries
import numpy
import Image

# Imports from brp plotting library
from brp.svg.base import SVGCanvas, PlotContainer, TextFragment
from brp.svg.plotters.line import LinePlotter
from brp.svg.plotters.raster import RasterPlotterMixin
from brp.svg.plotters.gradient import RGBGradient, GradientPlotter
from brp.svg.plotters.histogram import HistogramPlotter
from brp.svg.plotters.limit import YLimitPlotter

# Imports from the ssps single-pulse search library
from ssps import inf
from ssps.candidate import read_delays, SinglePulseReaderCondensed
from ssps.support import check_delays_option

# =============================================================================
# == Temporary copy of utility functions from other single pulse scripts,    ==
# == to be replaced with library functions when everything is public.        ==
# =============================================================================
SP_PATTERN = re.compile(r'\S+_DM(?P<dm>\d+\.\d+)\.singlepulse')
INF_PATTERN = re.compile(r'\S+_DM(?P<dm>\d+\.\d+)\.inf')


def find_files(directory, pattern):
    '''
    Find files matching pattern, return map from DM to filename.
    '''

    files = os.listdir(directory)
    dm2file = {}
    for f in files:
        m = pattern.match(f)
        if m:
            dm = float(m.group('dm'))
            dm2file[dm] = f

    return dm2file


def read_metadata(directory, dm2inf):
    '''
    Read all given PRESTO style .inf files in directory.
    '''

    metadata_map = {}
    for dm, f in dm2inf.iteritems():
        metadata_map[dm] = inf.inf_reader(os.path.join(directory, f))

    return metadata_map


# =============================================================================
# =============================================================================


def get_commandline_parser():
    '''
    Set up commandline parser.
    '''
    p = optparse.OptionParser()
    p.add_option('--dmspercell', type='int', default=5, metavar='NUMBER',
                 dest='dmspercell',
                 help='Extend in number of DM trials for each grid cell.')
    p.add_option('--o', type='string', default='condensed', metavar='BASENAME',
                 dest='basename', help='Basename for output files.')
    p.add_option('--odir', type='string', metavar='OUT_DIR',
                 dest='odir', help='Output directory.')
    p.add_option('--s', type='float', default=0, metavar='TIME',
                 dest='s',
                 help='Start time in seconds, time as is in .singlepulse files.')
    p.add_option('--e', type='float', default=3000, metavar='TIME',
                 dest='e',
                 help='End time in seconds, time as is in .singlepulse files.')
    p.add_option('--black', action='store_true', dest='black', default=False,
                 help='Use black background for main panel (off by default).')
    p.add_option('--overwrite', action='store_true', dest='overwrite',
                 default=False,
                 help='Overwrite existing plots (default False).')
    p.add_option('--lo', dest='lo', type='float', metavar='LO_DM',
                 help='Lowest DM to plot (default is lowest from data).')
    p.add_option('--hi', dest='hi', type='float', metavar='HI_DM',
                 help='Highest DM to plot (default is highest from data).')
    p.add_option('--secondspercell', dest='secondspercell', metavar='SECONDS',
                 help='Number of seconds per grid cell.', type=float,
                 default=10)
    p.add_option('--limit', dest='limit', type='float', metavar='N_DETECTIONS',
                 help='Draw line in lower panel at ... detections per bin.')
    p.add_option('--uselinkplaceholder', dest='uselinkplaceholder',
                 help='Use placeholders for next, previous and home links.',
                 default=False, action='store_true')
    p.add_option('--type', dest='type', type='string', metavar='TYPE_OF_PLOT',
                 help='Binning type: N number detections, M max(snr), S sum(snr)',
                 default='N')
    p.add_option('--delays', dest='delays', type='string',
                 help='File with DM delays (two columns; <DM> <delay(s)>).')

    return p


def check_type_option(options, args, p):
    '''
    Check --type commandline option.
    '''
    ALLOWED_TYPES = 'NMS'
    if options.type not in ALLOWED_TYPES:
        print 'Choose a type of plot from %s' % str(ALLOWED_TYPES.split())
        p.print_usage()
        p.print_help()
        sys.exit(1)


def get_dmi_range(spr, dmspercell):

    if len(spr.dms) % dmspercell == 0:
        max_dmi = len(spr.dms) - 1
    else:
        max_dmi = dmspercell * (1 + len(spr.dms) // dmspercell) - 1

    return 0, max_dmi


def count_detections(spr, dmspercell, max_dmi, start_time, end_time,
                     secondspercell, plot_type):
    ybins = (max_dmi + 1) // dmspercell
    assert (max_dmi + 1) % dmspercell == 0

    xbins = (end_time - start_time) / secondspercell
    if xbins != int(xbins):
        xbins = int(xbins) + 1
    dx = secondspercell

    ar = numpy.zeros((xbins, ybins), dtype=numpy.dtype(int))
    for i, dm in enumerate(spr.dms):
        ycell = i // dmspercell

        # TODO : see whether the delay compensation can be moved here?!
        if plot_type == 'N':
            for t, snr in spr.iterate_trial(dm):
                if start_time <= t < end_time:
                    xcell = int((t - start_time) / dx)
                    ar[xcell, ycell] += 1
        elif plot_type == 'M':
            for t, snr in spr.iterate_trial(dm):
                if start_time <= t < end_time:
                    xcell = int((t - start_time) / dx)
                    ar[xcell, ycell] = max(snr, ar[xcell, ycell])
        elif plot_type == 'S':
            for t, snr in spr.iterate_trial(dm):
                if start_time <= t < end_time:
                    xcell = int((t - start_time) / dx)
                    ar[xcell, ycell] += snr

    return ar


def colorcoded_ar_2d2png_string(color_ar):
    # Swap the axes of array so that the resulting PIL image comes out
    # oriented correctly.
    color_ar = numpy.swapaxes(color_ar, 0, 1)
    color_ar = numpy.flipud(color_ar)
    # Create a PNG image from the histogram
    im = Image.fromarray(color_ar)
    tmp = StringIO.StringIO()
    im.save(tmp, format='png')
    return 'data:image/png;base64,\n' + encodestring(tmp.getvalue())


def colorcode_ar_2d(ar, gradient):
    shape = ar.shape
    colors = numpy.zeros((shape[0], shape[1], 3), dtype=numpy.uint8)

    if numpy.amax(ar) != 0:
        for ix in range(shape[0]):
            for iy in range(shape[1]):
                r, g, b = gradient.get_color(ar[ix, iy])
                colors[ix, iy] = (r * 255, g * 255, b * 255)
    return colors


class SinglePulsePlotter(RasterPlotterMixin):
    def __init__(self, png_str, bbox):
        self.encoded_png = png_str
        self.img_bbox = bbox


def check_directories(options, args, p):
    output_files = []
    for plot_i, searchoutdir in enumerate(args):
        fn = options.basename + '-%08d.xml' % plot_i
        if options.odir is not None:
            print 'With outputdir %s' % options.odir
            if not os.path.exists(options.odir):
                tmp = os.path.realpath(options.odir)
                os.makedirs(tmp)
                fn = os.path.join(tmp, fn)
            else:
                fn = os.path.join(options.odir, fn)

        if os.path.exists(fn):
            if not options.overwrite:
                raise Exception('Won\'t overwrite %s' % fn)
        output_files.append((plot_i, searchoutdir, fn))
    return output_files


if __name__ == '__main__':
    print 'Called with:', sys.argv

    p = get_commandline_parser()
    options, args = p.parse_args()

    if len(args) == 0:
        p.print_usage()
        p.print_help()

    # Some custom checks:
    delays_file = check_delays_option(options, args, p)
    check_type_option(options, args, p)
    output_files = check_directories(options, args, p)

    for plot_i, searchoutdir, fn in output_files:
        print plot_i, searchoutdir, fn
        # Do some setting up and calculations:
        cv = SVGCanvas(1250, 760)
        try:
            print 'Scanning for .singlepulse and .inf files ...'
            spr = SinglePulseReaderCondensed(searchoutdir, options.s, options.e,
                                        delays_file, options.lo, options.hi)
            print '... done'
        except:
            raise
            datapath = os.path.abspath(searchoutdir)
            msg = 'Problem with data in %s, nothing present?' % datapath
            cv.add_plot_container(TextFragment(100, 100, msg, font_size=15))
        else:
            print 'Number of DM-trials: %d' % len(spr.dms)
            tf = TextFragment(870, 600, spr.sp_map[spr.dms[0]][:-12])
            cv.add_plot_container(tf)
            min_dmi, max_dmi = get_dmi_range(spr, options.dmspercell)

            # for main, count detections (or do SNR calculation per cell):
            ar = count_detections(spr, options.dmspercell, max_dmi, options.s,
                                  options.e, options.secondspercell,
                                  options.type)
            # set up the gradient:
            if options.type == 'N':
                m, M = 1, 30
            elif options.type == 'M':
                m, M = 5, 30
            elif options.type == 'S':
                m, M = 5, 30

            if options.black:
                gr = RGBGradient((m, M), (0, 0, 1), (1, 0, 0), min_value=1,
                                 max_value=40)
            else:
                gr = RGBGradient((m, M), (0, 0, 1), (1, 0, 0), min_value=1,
                                 max_value=40, min_value_color=(1, 1, 1),
                                 max_value_color=(1, 0, 0))

            colored_ar = colorcode_ar_2d(ar, gr)
            png_str = colorcoded_ar_2d2png_string(colored_ar)

            # for right panel
            if options.type in 'NS':
                collapsed_h = numpy.sum(ar, axis=0)
            elif options.type == 'M':
                collapsed_h = numpy.max(ar, axis=0)

            stepsize = options.dmspercell
            edges_dmi = [min_dmi + i * stepsize for i in
                         range(len(collapsed_h) + 1)]

            bins_dmi = []
            for ii, val in enumerate(collapsed_h):
                bins_dmi.append((edges_dmi[ii], edges_dmi[ii + 1], val))

            # bottom panel
            if options.type in 'NS':
                collapsed_v = numpy.sum(ar, axis=1)
            elif options.type == 'M':
                collapsed_v = numpy.max(ar, axis=1)

            dt = (options.e - options.s) / len(collapsed_v)
            edges_t = [options.s + i * dt for i in range(len(collapsed_v) + 1)]

            bins_t = []
            for ii, val in enumerate(collapsed_v):
                bins_t.append((edges_t[ii], edges_t[ii + 1], val))

            # Do the plots:
            PLOT_HEIGHT = 600

            TYPE2TITLE = {'N': 'Count', 'M': 'Max(SNR)', 'S': 'Sum(SNR)'}

            # Main panel, color-coded number of detections on DM-trial-time plane:
            pc_main = PlotContainer(0, -30, 900, PLOT_HEIGHT)
            pc_main.bottom.hide_label()
            pc_main.bottom.hide_tickmarklabels()
            pc_main.top.hide_label()
            pc_main.top.hide_tickmarklabels()
            pc_main.right.hide_tickmarklabels()
            pc_main.right.hide_label()
            pc_main.left.set_label('DM Index')
            pc_main.add_plotter(SinglePulsePlotter(png_str, (options.s, min_dmi,
                                options.e, max_dmi)))
            pc_main.set_minimum_data_bbox((options.s, min_dmi, options.e, max_dmi))
            cv.add_plot_container(pc_main)

            # First right-side panel, histogram of number of detections versus DM:
            pc_right1 = PlotContainer(820, -30, 270, PLOT_HEIGHT)
            pc_right1.add_plotter(HistogramPlotter(bins_dmi, False,
                                  orientation='vertical'))
            pc_right1.bottom.set_label(TYPE2TITLE[options.type])
            pc_right1.top.hide_label()
            pc_right1.top.hide_tickmarklabels()
            pc_right1.left.hide_tickmarklabels()
            pc_right1.left.hide_label()
            pc_right1.right.hide_tickmarklabels()
            pc_right1.right.hide_label()
            cv.add_plot_container(pc_right1)

            # Second right-side panel, DM-trial number versus DM:
            pc_right2 = PlotContainer(1010, -30, 270, PLOT_HEIGHT)
            pc_right2.left.hide_tickmarklabels()
            pc_right2.left.hide_label()
            pc_right2.right.hide_tickmarklabels()
            pc_right2.right.hide_label()
            pc_right2.top.hide_label()
            pc_right2.top.hide_tickmarklabels()
            pc_right2.bottom.set_label('DM')
            pc_right2.left.set_label('DM Index')
            dms = spr.dms[:]
            dm_indices = list(range(len(dms)))
            pc_right2.add_plotter(LinePlotter(dms, dm_indices,
                                  use_markers=False))
            pc_right2.set_minimum_data_bbox((dms[0], min_dmi, dms[1], max_dmi))
            cv.add_plot_container(pc_right2)

            # Bottom panel, histogram of number of detections versus t:
            pc_bottom = PlotContainer(0, PLOT_HEIGHT - 80 - 30, 900, 270)
            pc_bottom.add_plotter(HistogramPlotter(bins_t, False,
                                  orientation='horizontal'))
            if options.limit is not None:
                pc_bottom.add_plotter(YLimitPlotter(options.limit))
            pc_bottom.top.hide_tickmarklabels()
            pc_bottom.top.hide_label()
            pc_bottom.bottom.set_label('Time (s)')
            pc_bottom.right.hide_label()
            pc_bottom.left.set_label(TYPE2TITLE[options.type])
            pc_bottom.right.hide_tickmarklabels()
            pc_bottom.right.hide_label()
            cv.add_plot_container(pc_bottom)

            # Draw the gradient:
            grp = GradientPlotter(gr, 'horizontal')
            pc_gr = PlotContainer(820, PLOT_HEIGHT + 40, 460, 120, data_padding=0)
            pc_gr.left.hide_all()
            pc_gr.right.hide_all()
            pc_gr.top.hide_tickmarklabels()
            pc_gr.top.hide_label()
            pc_gr.bottom.set_label(TYPE2TITLE[options.type])
            pc_gr.add_plotter(grp)
            cv.add_plot_container(pc_gr)

        # Text (data set, next previous, etc...)
        if len(output_files) > 1:
            if plot_i < len(output_files) - 1:
                tf = TextFragment(870, 580, 'Next', color='blue',
                                  link=output_files[plot_i + 1][2])
                cv.add_plot_container(tf)
            if plot_i > 0:
                tf = TextFragment(920, 580, 'Previous', color='blue',
                                  link=output_files[plot_i - 1][2])
                cv.add_plot_container(tf)
        elif options.uselinkplaceholder:
                print 'Using link placeholders'
                cv.add_plot_container(TextFragment(870, 580, 'Home',
                                      color='blue',
                                      link='HOME_PLACEHOLDER'))
                cv.add_plot_container(TextFragment(1020, 580, 'Next',
                                      color='blue',
                                      link='NEXT_PLACEHOLDER'))
                cv.add_plot_container(TextFragment(1070, 580, 'Previous',
                                      color='blue',
                                      link='PREVIOUS_PLACEHOLDER'))

        with open(fn, 'w') as f:
            cv.draw(f)
        print spr.n_error, spr.n_success
