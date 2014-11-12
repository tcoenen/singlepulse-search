#!/usr/bin/env python
'''
Find pulses, sort by SNR, look for double detections.

Note: built using the new single pulse grouping algorithm. This script will
in future be filtering out bad data on the fly, allowing it to be run on bad
messy data sets.
'''
import os
import sys
import optparse
import time

import ssps
from ssps import candidate
print candidate.__file__
from ssps.candidate.vanilla import SinglePulseReaderBase
from ssps import pulse
from ssps import pulsetrain
from ssps import diagnostic
from ssps.support import check_delays_option

# -----------------------------------------------------------------------------
# -- Settings for data reduction pipeline -------------------------------------

DMS_ADJACENT = 4

# Commandline help message:
USAGE = '''python grab.py <options>

Grab promising single pulses from relatively clean single pulse data sets and
create diagnostic plots showing them in red on a background of black single
pulse detections.

In the output directory the following files will be present:
    <basename>_DM<dm>_<number>.xml: An SVG file with the plot mentioned above.
    <basename>_DM<dm>_<number>.peaks.singlepulsee: A text file with the DM,
        time, SNR of the interesting pulses.
'''


def get_commandline_parser():
    '''
    Set up commandline parser.
    '''
    p = optparse.OptionParser(usage=USAGE)
    p.add_option('--searchoutdir', dest='searchoutdir', type='string',
                 help='Search output dir (with INF and SINGLEPULSE subdir).')
    p.add_option('--odir', dest='odir', type='string',
                 help='Output directory for results of search.')
    p.add_option('--o', dest='o', type='string',
                 help='Basename of output files')
    p.add_option('--snr', dest='snr', type='float', default=7,
                 help='Minimum peak signal-to-noise for bright single pulse.')
    p.add_option('--repeats', dest='min_n_repeats', type='int', default=5,
                 help='Minimum no. of bright single pulses at a certain DM.')
    p.add_option('--tlo', dest='tlo_dm', type='float',
                 help='Minimum trigger DM to search (default DM=0).')
    p.add_option('--thi', dest='thi_dm', type='float',
                 help='Maximum trigger DM to search (no max by default)')
    p.add_option('--lo', dest='lo_dm', type='float', default=0,
                 help='Minimum DM to read in (default DM=0).')
    p.add_option('--hi', dest='hi_dm', type='float',
                 help='Maximum DM to read in (no max by default)')
    p.add_option('--ndmsbright', dest='n_dms_bright', type='int', default=6,
                 help='No. of DMs to search near bright pulse (default 6).')
    p.add_option('--ncandbright', dest='min_n_candidates_bright', type='int',
                 default=7,
                 help='Minimum no. of candidates per bright pulse (default 7)')
    p.add_option('--ncanddim', dest='min_n_candidates_dim', type='int',
                 default=5,
                 help='Minimum no. of candidates per dim pulse (deafult 5).')
    p.add_option('--ndms', dest='ndms', type='int',
                 help='Minimum no. of DMs present in data set (sanity check).')
    p.add_option('--overwrite', action='store_true', dest='overwrite',
                 default=False,
                 help='Overwrite existing plots (default False).')
    p.add_option('--uselinkplaceholder', dest='uselinkplaceholder',
                 help='Use placeholders for next, previous and home links.',
                 default=False, action='store_true')
    p.add_option('--marker', dest='marker', type='string',
                 help='Mark DMs (comma separated list, no spaces).')
    p.add_option('--delays', dest='delays', type='string',
                 help='File with DM delays (two columns; <DM> <delay(s)>).')

    return p


def check_marker_option(options, args, p):
    '''
    Check --marker commandline option.
    '''
    marker_dms = []
    if options.marker is not None:
        chunks = options.marker.split(',')
        for chunk in chunks:
            if not chunk:
                continue
            try:
                dm = float(chunk)
            except:
                raise Exception('--marker incorrectly specified')
            else:
                marker_dms.append(dm)

    return marker_dms


def check_searchoutdir_option(options, args, p):
    '''
    Check --searchoutdir commandline option.
    '''
    if  options.searchoutdir is None:
        print 'No search output directory provided use --searchoutdir option!'
        p.print_help()
        sys.exit(1)
    else:
        searchoutdir = os.path.abspath(options.searchoutdir)

    return searchoutdir


def check_odir_option(options, args, p):
    '''
    Check --odir commandline option.
    '''
    if options.odir is None:
        print 'No output directory provided for results use --odir option!'
        p.print_help()
        sys.exit(1)
    else:
        OUTDIR = os.path.abspath(options.odir)
        if not os.path.exists(OUTDIR):
            print 'Output directory does not exist! (Exiting!)'
            p.print_help()
            sys.exit(1)

    return OUTDIR


def handle_commandline():
    '''
    Set up commandline parser, check and return options.
    '''

    p = get_commandline_parser()
    options, args = p.parse_args()

    # custom checks:
    marker_dms = check_marker_option(options, args, p)
    delays_file = check_delays_option(options, args, p)
    searchoutdir = check_searchoutdir_option(options, args, p)
    odir = check_odir_option(options, args, p)

    return options, args, marker_dms, delays_file, searchoutdir, odir

if __name__ == '__main__':
    # Echo the commandline for crude logging:
    print 'Called with :', sys.argv
    t0 = time.time()
    # Handle the scripts commandline options:
    options, args, marker_dms, delays_file, searchoutdir, odir = \
        handle_commandline()
    # read data
    print '=' * 77
    print 'Processing %s' % searchoutdir
    print 'Looking for datafiles.'
    spr = SinglePulseReaderBase(
        searchoutdir, 
        None, # tstart
        None, # tend
        delays_file, 
        lodm=options.lo_dm,
        hidm=options.hi_dm,
        max_downfact=30
    )

    print 'DM range after reading: [%.2f, %.2f]' % (spr.dms[0], spr.dms[-1])

    if options.ndms is not None:
        ndms = len(spr.dms)
        if ndms < options.ndms:
            print 'Too few DM trials, aborting!'
            print 'Needed %d and found %d DM trials.' % (options.ndms, ndms)
            sys.exit(1)
    # determine the basename for the output files:
    if options.o is None:
        dm_str = '_DM%.2f' % spr.dms[0]
        basename = spr.md_map[spr.dms[0]].data_file[:-len(dm_str)]
    else:
        basename = options.o
    print 'Looking for pulses.'
    # Call the rewritten ssps candidate grouping algorithm to find single
    # pulses by combining candidates across DM trials.
    pulses = pulse.group(spr, DMS_ADJACENT)
    pulses = pulse.annotate_pulses(spr, pulses)

    # Find the bright pulses and the dim ones, sort them for SNR.
    # find pulsar/RRAT candidates
    print 'Sifting groups, single pulse detections, into bright and dim ones.'
    print '  Threshold (peak) snr for bright pulses %f .' % options.snr
    print '  Minimum number of candidates per group %d .' % \
        options.min_n_candidates_bright
    print '  Minimum no. of bright pulses to consider a DM interesting %d.' % \
        options.min_n_repeats
    pulse_trains, rejects = \
        pulsetrain.extract_pulsetrains(spr, pulses, options.snr,
                                       options.min_n_repeats,
                                       options.min_n_candidates_bright,
                                       options.n_dms_bright,
                                       options.tlo_dm, options.thi_dm)

    n_bright = sum(len(x) for x in pulse_trains)
    print '  Found %d bright pulses near %d trial DMs.' % \
        (n_bright, len(pulse_trains))
    print 'Adding dim pulses to the bright ones found already.'
    print '  Minimum number of candidates per group %d .' % \
        options.min_n_candidates_dim
    pulse_trains, rejects = \
        pulsetrain.add_dim_pulses(spr, pulse_trains, rejects,
                                  options.min_n_candidates_dim)
    n_dim = sum(len(x) for x in pulse_trains) - n_bright
    print '  Found %d dim pulses.' % n_dim

    for keepers in pulse_trains:
        print 'DM %.2f had %d detections with max(snr) = %.2f.' % \
            (keepers[0].dm, len(keepers), max(x.snr for x in keepers))

    print '''TOOK %.2f SECONDS.''' % (time.time() - t0)

    # Determine the output file names:
    basenames = []
    for i, keepers in enumerate(pulse_trains):
        fn = basename + ('_%04d' % i) + ('_DM%.2f' % (keepers[0].dm))
        if not options.overwrite:
            if os.path.exists(os.path.join(options.odir, fn + '.xml')) or \
                    os.path.exists(os.path.join(options.odir,
                                   fn + '.peaks.singlepulse')):
                print '\nERROR: Output files exist, refusing overwrite!'
                sys.exit(1)

        basenames.append(fn)

    # Plot everything:
    n_plots = len(pulse_trains)
    for i in range(n_plots):
        # Determine the full path to this diagnostic plot:
        fn = os.path.join(options.odir, basenames[i] + '.xml')
        print 'Writing plot to %s' % fn
        # Determine what the links in the diagnostic plot should point to:
        if options.uselinkplaceholder:
            PREVIOUS_LINK = 'PREVIOUS_PLACEHOLDER'
            NEXT_LINK = 'NEXT_PLACEHOLDER'
        else:  # Remember, we want relative links here!
            if i > 0:
                PREVIOUS_LINK = basenames[i - 1] + '.xml'
            else:
                PREVIOUS_LINK = ''
            if i < n_plots - 1:
                NEXT_LINK = basenames[i + 1] + '.xml'
            else:
                NEXT_LINK = ''
        # Make the plot:
        diagnostic.plot(fn, spr, options, pulse_trains, i,
                        rejects, NEXT_LINK, PREVIOUS_LINK, marker_dms)
        # Write a file with the peaks of each pulse in the pulsetrain (in the
        # same format as the PRESTO .singlepulse files).
        pulse.write_arrivals(os.path.join(options.odir,
                             basenames[i] + '.peaks.singlepulse'), pulse_trains[i])
    print 'Took %.2f seconds!' % (time.time() - t0)
