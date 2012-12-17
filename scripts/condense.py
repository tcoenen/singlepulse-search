#!/usr/bin/env python
# By Thijs Coenen (2012) for PhD thesis research at the Universiteit van 
# Amsterdam.
'''
Condensed single pulse plots from PRESTO .singlepulse/.inf files.

Second large iteration, meant to produce plots that can be included in thesis,
so no black backgrounds and glitches fixed.
'''
# Python standard library imports
from __future__ import division
import os
import optparse
import StringIO
from base64 import encodestring
import re

# Standard 3rd party libraries
import numpy
import Image

# Imports from brp plotting library
from brp.svg.base import SVGCanvas, PlotContainer, TextFragment
from brp.svg.plotters.scatter import ScatterPlotter
from brp.svg.plotters.raster import RasterPlotterMixin
from brp.svg.plotters.gradient import RGBGradient, GradientPlotter
from brp.svg.plotters.histogram import HistogramPlotter
from brp.svg.plotters.limit import YLimitPlotter

# Other single pulse search imports
import inf

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


def check_commandline():
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
    options, args = p.parse_args()
    return options, args


class SinglePulseReaderBase(object):
    def __init__(self, directory, tstart, tend, lodm=None, hidm=None):
        self.sp_dir = os.path.join(directory, 'SINGLEPULSE')
        self.inf_dir = os.path.join(directory, 'INF')
        self.sp_map = find_files(self.sp_dir, SP_PATTERN)
        self.inf_map = find_files(self.inf_dir, INF_PATTERN)
        self.md_map = read_metadata(self.inf_dir, self.inf_map)
        self.dms = list(set(self.sp_map.keys()) & set(self.inf_map.keys()))

        if len(self.dms) == 0:
            raise Exception('No matching .inf AND .singlepulse files in %s' %
                            directory)

        # Store desired DM, time range.
        self.tstart = tstart
        self.tend = tend
        self.lodm = lodm
        self.hidm = hidm

        # Only work on the DM trials in the desired range -> filter them.
        if self.lodm:
            self.dms = [dm for dm in self.dms if lodm <= dm]
        if self.hidm:
            self.dms = [dm for dm in self.dms if dm <= hidm]

        self.dms.sort()
        self.dm2idx = dict([(dm, i) for i, dm in enumerate(self.dms)])

        self.n_success = 0
        self.n_error = 0
        self.n_rejected = 0

    def iterate_trial(self, dm):
        # for each candidate return: t
        sp_file = os.path.join(self.sp_dir, self.sp_map[dm])

        with open(sp_file, 'r') as f:
            for line in f:
                split_line = line.split()
                try:
                    t = float(split_line[2])
                except:
                    self.n_error += 1
                else:
                    self.n_success += 1
                    if self.tstart <= t <= self.tend:
                        yield t


def get_dmi_range(spr, dmspercell):
    # Remember zero based indexing, hence the - 1 !
    if len(spr.dms) % dmspercell == 0:
        max_dmi = len(spr.dms) - 1
    else:
        max_dmi = dmspercell * (1 + len(spr.dms) // dmspercell) - 1

    return 0, max_dmi


def count_detections(spr, dmspercell, max_dmi, start_time, end_time,
                     secondspercell):
    ybins = (max_dmi + 1) // dmspercell
    assert (max_dmi + 1) % dmspercell == 0

    xbins = (end_time - start_time) / secondspercell
    if xbins != int(xbins):
        xbins = int(xbins) + 1
    dx = secondspercell

    ar = numpy.zeros((xbins, ybins), dtype=numpy.dtype(int))
    for i, dm in enumerate(spr.dms):
        ycell = i // dmspercell
        for t in spr.iterate_trial(dm):
            if start_time <= t < end_time:
                xcell = int((t - start_time) / dx)
                ar[xcell, ycell] += 1

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


def check_directories(dirs):
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
    options, args = check_commandline()

    output_files = check_directories(args)
    for plot_i, searchoutdir, fn in output_files:
        print plot_i, searchoutdir, fn
        # Do some setting up and calculations:
        cv = SVGCanvas(1250, 760)
        try:
            spr = SinglePulseReaderBase(searchoutdir, options.s, options.e,
                                        options.lo, options.hi)
        except:
            datapath = os.path.abspath(searchoutdir)
            msg = 'Problem with data in %s, nothing present?' % datapath
            cv.add_plot_container(TextFragment(100, 100, msg, font_size=15))
        else:
            print 'Number of DM-trials: %d' % len(spr.dms)
            tf = TextFragment(870, 580, spr.sp_map[spr.dms[0]][:-12])
            cv.add_plot_container(tf)
            min_dmi, max_dmi = get_dmi_range(spr, options.dmspercell)

            # for main panel:
            ar = count_detections(spr, options.dmspercell, max_dmi, options.s,
                                  options.e, options.secondspercell)
            if options.black:
                gr = RGBGradient((1, 30), (0, 0, 1), (1, 0, 0), min_value=1,
                                 max_value=40)
            else:
                gr = RGBGradient((1, 30), (0, 0, 1), (1, 0, 0), min_value=1,
                                 max_value=40, min_value_color=(1, 1, 1),
                                 max_value_color=(1, 0, 0))
            colored_ar = colorcode_ar_2d(ar, gr)
            png_str = colorcoded_ar_2d2png_string(colored_ar)

            # for right panel
            collapsed_h = numpy.sum(ar, axis=0)
            stepsize = options.dmspercell
            edges_dmi = [min_dmi + i * stepsize for i in
                         range(len(collapsed_h) + 1)]

            bins_dmi = []
            for ii, val in enumerate(collapsed_h):
                bins_dmi.append((edges_dmi[ii], edges_dmi[ii + 1], val))

            # bottom panel
            collapsed_v = numpy.sum(ar, axis=1)
            dt = (options.e - options.s) / len(collapsed_v)
            edges_t = [options.s + i * dt for i in range(len(collapsed_v) + 1)]

            bins_t = []
            for ii, val in enumerate(collapsed_v):
                bins_t.append((edges_t[ii], edges_t[ii + 1], val))

            # Do the plots:
            PLOT_HEIGHT = 600

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
            pc_right1.bottom.set_label('Count')
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
            pc_right2.add_plotter(ScatterPlotter(dms, dm_indices))
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
            pc_bottom.left.set_label('Count')
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
            pc_gr.bottom.set_label('Count')
            pc_gr.add_plotter(grp)
            cv.add_plot_container(pc_gr)

        # Text (data set, next previous, etc...)
        if len(output_files) > 0:
            if plot_i < len(output_files) - 1:
                tf = TextFragment(870, 600, 'Next', color='blue',
                                  link=output_files[plot_i + 1][2])
                cv.add_plot_container(tf)
            if plot_i > 0:
                tf = TextFragment(920, 600, 'Previous', color='blue',
                                  link=output_files[plot_i - 1][2])
                cv.add_plot_container(tf)

        with open(fn, 'w') as f:
            cv.draw(f)
