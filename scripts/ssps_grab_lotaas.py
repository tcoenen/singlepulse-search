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

#from ssps import candidate
from ssps import pulse
from ssps import pulsetrain
from ssps import diagnostic
from ssps.support import check_delays_option

# for the porting:
from ssps import inf
import copy
from collections import defaultdict

# -----------------------------------------------------------------------------
# -- Settings for data reduction pipeline -------------------------------------

HARD_LIMIT_N_CANDIDATES = 1000000  # Maximum number of candidates in data.
EPSILON = 1e-5  # (Smaller than DM0.00 time-series bin size in seconds!)


class SinglePulseReaderLOTAAS(object):
    def __init__(self, basename, tstart, tend, delays_file, lodm=None, hidm=None,
                 max_downfact=30):
        # files to read:
        self.singlepulse_file = basename + '.singlepulse'
        self.inf_file = basename + '.inf'

        index, binwidth_map = self.scan_data()
        self.index = index
        self.dms = index.keys()
        self.dms.sort()

        metadata = inf.inf_reader(self.inf_file)
        metadata_map = {}
        for dm in self.dms:
            metadata_map[dm] = copy.deepcopy(metadata)
            # LOTAAS specific HACK to get around missing .inf files: 
            metadata_map[dm].binwidth = binwidth_map[dm] 
        self.md_map = metadata_map

        delay_map = defaultdict(float)
        if delays_file:
            print 'Loading delays from %s' % delays_file
            delay_map = read_delays(delays_file, delay_map)
        self.delay_map = delay_map

        self.dm2idx = dict([(dm, i) for i, dm in enumerate(self.dms)])
        self.max_downfact = max_downfact

        self.n_success = 0
        self.n_error = 0
        self.n_rejected = 0
        
        self.tstart = tstart
        self.tend = tend

        # Only work on the DM trials in the desired range -> filter them.
        if lodm:
            self.dms = [dm for dm in self.dms if lodm <= dm]
        if hidm:
            self.dms = [dm for dm in self.dms if dm <= hidm]

        self.dm_delay_map = dict((dm, 0) for dm in self.dms)

    def scan_data(self):
        '''
        Extract binwidths and DM positions from the singlepulse file.

        Note:
        Assumes that the LOTAAS style singlepulse file contains many DMs and that
        the data for each DM is a set of consecutive lines that are sorted in time.

        '''
        index = {}
        binwidth_map = {}
        last_dm = -1

        with open(self.singlepulse_file, 'r') as f:
            filepos = 0
            line = f.readline()

            while len(line) > 0:
                lastpos = filepos
                filepos = f.tell()

                if line and line[0] != '#':
                    split_line = line.split()
                    dm = float(split_line[0])
                    
                    if dm != last_dm:
                        if last_dm != -1:
                            index[last_dm][1] = lastpos
                        index[dm] = [lastpos, filepos]
                        binwidth_map[dm] = float(split_line[5])

                    last_dm = dm
                else:
                    assert lastpos == 0  # only want comments on first line of file!

                line = f.readline()

        return index, binwidth_map

    def get_t_overlap(self, dm):
        '''
        TBD
        '''
        return self.max_downfact * self.md_map[dm].bin_width + EPSILON

    def iterate_trial(self, dm):

        epsilon = EPSILON
        max_n_candidates = HARD_LIMIT_N_CANDIDATES
        bin_width = self.md_map[dm].bin_width
        i = self.dm2idx[dm]
        dm_delay = self.dm_delay_map[dm]

        delay = self.delay_map[dm]
        startpos, endpos = self.index[dm]

        with open(self.singlepulse_file, 'r') as f:
            f.seek(startpos)
            
            while f.tell() < endpos:
                line = f.readline()
                split_line = line.split()

                try:
                    snr = float(split_line[1])
                    t = float(split_line[2]) + dm_delay
                    sample = int(split_line[3])
                    downfact = int(split_line[4])
                    delta_t = downfact * bin_width + epsilon
                except:
                    self.n_error += 1
                else:
#                    if self.is_ok(t, dm):
                    if True:
                        cand = (dm, snr, t, sample, downfact, t - delta_t,
                                t + delta_t, i)

                        self.n_success += 1
                        if self.n_success > max_n_candidates:
                            raise TooManyCandidates(max_n_candidates)
                        yield cand
                    else:
                        self.n_rejected += 1


def get_dmi_range(spr, dmspercell):

    if len(spr.dms) % dmspercell == 0:
        max_dmi = len(spr.dms) - 1
    else:
        max_dmi = dmspercell * (1 + len(spr.dms) // dmspercell) - 1

    return 0, max_dmi

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


DMS_ADJACENT = 4

# Commandline help message:
USAGE = '''python grab.py <options>

Grab promising single pulses from relatively clean single pulse data sets and
create diagnostic plots showing them in red on a background of black single


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
    odir = check_odir_option(options, args, p)

    return options, args, marker_dms, delays_file, odir

if __name__ == '__main__':
    # Echo the commandline for crude logging:
    print 'Called with :', sys.argv
    t0 = time.time()
    # Handle the scripts commandline options:
    options, args, marker_dms, delays_file, odir = \
        handle_commandline()
    # read data
    print '=' * 77
#    print 'Processing %s' % searchoutdir
    print 'Looking for datafiles.'
    spr = SinglePulseReaderLOTAAS(args[0], delays_file, 30, options.lo_dm,
                                  options.hi_dm)

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
