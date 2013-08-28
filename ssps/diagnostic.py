#!/usr/bin/env python
from __future__ import division

from brp.svg.base import SVGCanvas, PlotContainer, TextFragment
from brp.svg.plotters.histogram import HistogramPlotter
from brp.svg.plotters.scatter import ScatterPlotter
from brp.svg.plotters.limit import YLimitPlotter
from inf import dec2ascii, ra2ascii

from brp_diagnostic import DetectionPlotterTuples

PADDING = 10
MAX_N_CANDIDATES = 20000  # Fall back to rasterized beyond this number of cands

# -----------------------------------------------------------------------------
# -- General helper functions -------------------------------------------------


def find_dm_range(dms, pulsetrain):
    '''Find DM range of candidates in pulsetrain.'''
    # TODO : see whether this works better if we work in DM index space
    # in stead of DM space.
    lo_dm = pulsetrain[0].body[0][0]
    hi_dm = pulsetrain[0].body[0][0]

    for pulse in pulsetrain:
        if pulse.lo_dm < lo_dm:
            lo_dm = pulse.lo_dm
        if pulse.hi_dm > hi_dm:
            hi_dm = pulse.hi_dm

    # Make the range of DMs we plot 3 times a wide (if possible):
    w_dm = hi_dm - lo_dm
    mean_dm = (hi_dm + lo_dm) / 2
    hi_dm = min(mean_dm + w_dm, dms[-1])
    lo_dm = max(dms[0], mean_dm - w_dm)

    return lo_dm, hi_dm


def flatten_pulses(pulses, lo_dm, hi_dm):
    '''Make a flat list of single pulse candidates to be plot.'''
    candidates = []
    for pulse in pulses:
        candidates.extend(c for c in pulse.body if lo_dm <= c[0] <= hi_dm)

    return candidates

# -----------------------------------------------------------------------------
# -- Main panel plotting function ---------------------------------------------


def plot_main_panel(cv, foreground, background, lodm, hidm, marker_dms,
                    fallback):
    '''
    Add main panel to pulsetrain diagnostic plot.

    Note: for internal use, call the diagnostic.plot() function in stead.
    '''
    # Main panel, detections on the time-DM plane
    pc1 = PlotContainer(0, 150, 830, 600, data_padding=PADDING)
    pc1.right.hide_label()
    pc1.right.hide_tickmarklabels()
    pc1.bottom.hide_label()
    pc1.left.set_label('DM')
    pc1.top.set_label('Time (s)')
    pc1.add(DetectionPlotterTuples(background, color='black'), fallback)
    pc1.add(DetectionPlotterTuples(foreground, color='red'), fallback)
    pc1.set_minimum_y_range(lodm, hidm)

    for m_dm in marker_dms:
        if lodm <= m_dm <= hidm:
            pc1.add(YLimitPlotter(m_dm, color='orange'))

    cv.add(pc1)

# -----------------------------------------------------------------------------
# -- First right-hand side panel ----------------------------------------------


def plot_snr_scatter(cv, foreground, background, lodm, hidm, marker_dms,
                     fallback):
    '''
    Add SNR versus DM scatter plot pulsetrain diagnostic plot.

    Note: for internal use, call the diagnostic.plot() function in stead.
    '''
    pc2 = PlotContainer(740, 150, 300, 600, data_padding=PADDING)
    pc2.right.hide_label()
    pc2.right.hide_tickmarklabels()
    pc2.left.hide_label()
    pc2.left.hide_tickmarklabels()
    pc2.bottom.hide_label()
    pc2.top.set_label('SNR')
    scp2 = ScatterPlotter([x[1] for x in foreground],
                          [x[0] for x in foreground], color='red')
    scp3 = ScatterPlotter([x[1] for x in background],
                          [x[0] for x in background], color='black')
    pc2.add(scp3, fallback)
    pc2.add(scp2, fallback)
    pc2.set_minimum_y_range(lodm, hidm)
    for m_dm in marker_dms:
        if lodm <= m_dm <= hidm:
            pc2.add(YLimitPlotter(m_dm, color='orange'))
    cv.add(pc2)

# -----------------------------------------------------------------------------
# -- Second right-hand side panel ---------------------------------------------


def count_candidates(dms, candidates):
    '''
    Count single pulse candidates on grid, return dm -> count dictionary.
    '''
    count_dict = dict((dm, 0) for dm in dms)
    for c in candidates:
        try:
            count_dict[c[0]] += 1  # read [0] -> .dm
        except KeyError, e:
            print c[0], c
            raise e
    counts = [count_dict[dm] for dm in dms]
    return counts


def bin_candidates(all_dms, foreground_cands, background_cands, lo_dm, hi_dm):
    '''Create number of candidates as function of DM histograms.'''
    # For reference, plotting library 'brp' uses the following layout for
    # histograms: [(x1, x2, value), ...] (Which is what we construct here.)

    # Determine the bins we need (aligned with DM trials):
    dms = [dm for dm in all_dms if lo_dm <= dm <= hi_dm]
    midpoints = [(dms[i] + dms[i + 1]) / 2 for i in range(len(dms) - 1)]
    edges = [dms[0]]
    edges.extend(midpoints)
    edges.append(dms[-1])

    # Count the foreground, background and total:
    f_counts = count_candidates(dms, foreground_cands)
    b_counts = count_candidates(dms, background_cands)

    n_dms = len(dms)
    f_bins = [(edges[i], edges[i + 1], f_counts[i]) for i in range(n_dms)]
    b_bins = [(edges[i], edges[i + 1], b_counts[i]) for i in range(n_dms)]
    t_bins = [(edges[i], edges[i + 1], b_counts[i] + f_counts[i])
              for i in range(n_dms)]

    return f_bins, b_bins, t_bins


def plot_count_histogram(cv, dms, foreground, background, lodm, hidm,
                         marker_dms):
    '''
    Add candidates per DM histogram to pulse train diagnostic.

    Note: for internal use, call the diagnostic.plot() function in stead.
    '''
    # Determine the brp style histogram bins:
    fg_bins, bg_bins, total_bins = \
        bin_candidates(dms, foreground, background, lodm, hidm)

    # Plot histogram of number of detections per DM trial
    pc3 = PlotContainer(950, 150, 300, 600, data_padding=PADDING)
    pc3.left.hide_label()
    pc3.left.hide_tickmarklabels()
    pc3.bottom.hide_label()
    pc3.top.set_label('N')
    pc3.add(HistogramPlotter(bg_bins, orientation='vertical', color='black'))
    pc3.add(HistogramPlotter(total_bins, orientation='vertical', color='gray'))
    pc3.add(HistogramPlotter(fg_bins, orientation='vertical', color='red'))
    pc3.right.set_label('DM')
    pc3.set_minimum_y_range(lodm, hidm)

    for m_dm in marker_dms:
        if lodm <= m_dm <= hidm:
            pc3.add(YLimitPlotter(m_dm, color='orange'))

    cv.add(pc3)

# -----------------------------------------------------------------------------
# -- Text for diagnostic plots ------------------------------------------------


def add_text(cv, metadata, next_link, previous_link, n_candidates, fallback):
    '''
    Add plot metadata and links.

    Note: for internal use, call the diagnostic.plot() function in stead.
    '''

    dm_str = '_DM%.2f' % metadata.dm
    basename = metadata.data_file[:-len(dm_str)]

    cv.add(TextFragment(50, 20, 'Right Ascension  %s' %
           ra2ascii(*metadata.j2000_ra)))
    cv.add(TextFragment(50, 40, 'Declination %s' %
           dec2ascii(*metadata.j2000_dec)))
    cv.add(TextFragment(50, 60, 'Epoch (MJD) %.15f' %
           metadata.epoch_mjd))
    if fallback:
        mode = 'raster'
    else:
        mode = 'vector'

    # links hack
    if next_link:
        cv.add(TextFragment(50, 100, 'Next', link=next_link, color='blue'))
    if previous_link:
        cv.add(TextFragment(50, 120, 'Previous', link=previous_link,
               color='blue'))
    # the data set
    cv.add(TextFragment(50, 140,
           '%s [%d candidates, %s mode]' % (basename, n_candidates, mode)))


def add_settings_text(cv, opts):
    '''
    Output settings used for pulsetrain extraction algorithm to plot.

    Note: for internal use, call the diagnostic.plot() function in stead.
    '''
    cv.add(TextFragment(450, 20,
           'Bright pulse peak SNR threshold: %.2f' % opts.snr))
    cv.add(TextFragment(450, 40,
           'Minimum bright pulses per pulse train: %d' % opts.min_n_repeats))
    cv.add(TextFragment(450, 60,
           'Minimum candidates per bright pulse: %d' %
           opts.min_n_candidates_bright))
    cv.add(TextFragment(450, 80,
           'Minimum candidates per dim pulse: %d' %
           opts.min_n_candidates_dim))
    cv.add(TextFragment(450, 100,
           'DM trials to search near bright pulse peak: %d' %
           opts.n_dms_bright))


def add_trigger_text(cv, pulse_trains, pt_idx):
    '''
    Add text to pulsetrain diagnostic giving trigger SNR and DM.

    Note: for internal use, call the diagnostic.plot() function in stead.
    '''
    trigger_snr = pulse_trains[pt_idx][0].snr
    trigger_dm = pulse_trains[pt_idx][0].dm

    tf2 = TextFragment(810, 220, 'SNR %.2f' % trigger_snr, color='red')
    tf3 = TextFragment(810, 240, 'DM %.2f' % trigger_dm, color='red')
    cv.add(tf2)
    cv.add(tf3)

# -----------------------------------------------------------------------------
# -- Plotting function that creates the overall diagnostic plot ---------------


def plot(filename, spr, options, pulse_trains, pt_idx, rejected_pulses,
         nlink, plink, marker_dms):
    '''
    Create a pulsetrain diagnostic plot.

    Note: Call this function to create the plot, not others in this module.
    '''
    # Find DM range for this plot:
    lodm, hidm = find_dm_range(spr.dms, pulse_trains[pt_idx])

    # Flatten pulses into list of candidates (foreground, background):
    foreground = flatten_pulses(pulse_trains[pt_idx], lodm, hidm)

    background = flatten_pulses(rejected_pulses, lodm, hidm)
    for i in range(len(pulse_trains)):
        if i != pt_idx:
            background.extend(flatten_pulses(pulse_trains[i], lodm, hidm))

    # Determine whether it is wise to switch to rasterized plotting (called
    # fallback).
    n_candidates = len(foreground) + len(background)
    if n_candidates > MAX_N_CANDIDATES:
        fallback = True
    else:
        fallback = False

    # Set up the overall canvas to wich all panels are plot:
    cv = SVGCanvas(1250, 750)
    # Add the three panels to the plot:
    plot_main_panel(cv, foreground, background, lodm, hidm, marker_dms,
                    fallback)
    plot_snr_scatter(cv, foreground, background, lodm, hidm, marker_dms,
                     fallback)
    plot_count_histogram(cv, spr.dms, foreground, background, lodm, hidm,
                         marker_dms)
    # Add the text identifying the data set etc. to the diagnostic plot:
    add_text(cv, spr.md_map[spr.dms[0]], nlink, plink, n_candidates, fallback)
    add_settings_text(cv, options)
    add_trigger_text(cv, pulse_trains, pt_idx)

    # Finally write everything out to an SVG (xml) file.
    with open(filename, 'w') as f:
        cv.draw(f)
