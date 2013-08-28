'''
Module that deals with PRESTO single pulse search output.

Contains:
    - utility functions that search for .singlepulse and .inf files
    - object that represents a single pulse data set
    - support for adding delays to individual trial-DMs

'''
from __future__ import division

import os
import re

import inf

# --------------------------------------------------------------------------
# -- Constants -------------------------------------------------------------

HARD_LIMIT_N_CANDIDATES = 1000000  # Maximum number of candidates in data.
EPSILON = 1e-5  # (Smaller than DM0.00 time-series bin size in seconds!)

# --------------------------------------------------------------------------
# -- Deal with directories full of single pulse detections -----------------

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

# --------------------------------------------------------------------------
# -- Support for extra delays ----------------------------------------------


class DelaysFileWrong(Exception):
    pass


def read_delays(filename, dm_delay_map):
    '''
    Read text file containing extra delays.

    Expected format for each line:
    <dm value> <delay in seconds>

    Lines starting with # are ignored.
    '''
    with open(filename, 'r') as f:
        for line in f:
            if line and line[0] != '#':
                split_line = line.split()
                try:
                    dm = float(split_line[0])
                    delay = float(split_line[1])
                except (ValueError, IndexError):
                    raise DelaysFileWrong('The delays file %s is wrong!'
                                          % filename)
                else:
                    dm_delay_map[dm] = delay
    return dm_delay_map

# --------------------------------------------------------------------------
# -- Read raw candidates ---------------------------------------------------


class TooManyCandidates(Exception):
    def __init__(self, n):
        self.msg = 'Data set contains too many candidates max = %d !' % n

    def __str__(self):
        return self.msg


class SinglePulseReaderBase(object):
    def __init__(self, directory, delays_file, max_downfact=30, lodm=None,
                 hidm=None):
        '''
        Read PRESTO single pulse data set.

        Note: current implementation assumes the .inf files are in a
        subdirectory of <directory> called INF and the .singlepulse files
        in a subdirectory of <directory> callsed SINGELPULSE.
        '''
        self.sp_dir = os.path.join(directory, 'SINGLEPULSE')
        self.inf_dir = os.path.join(directory, 'INF')
        self.sp_map = find_files(self.sp_dir, SP_PATTERN)
        self.inf_map = find_files(self.inf_dir, INF_PATTERN)
        self.md_map = read_metadata(self.inf_dir, self.inf_map)
        self.dms = list(set(self.sp_map.keys()) & set(self.inf_map.keys()))

        if len(self.dms) == 0:
            raise Exception('No matching .inf AND .singlepulse files in %s' %
                            directory)
        if lodm is not None:
            self.dms = [dm for dm in self.dms if dm >= lodm]
        if hidm is not None:
            self.dms = [dm for dm in self.dms if dm <= hidm]
        self.dms.sort()

        if len(self.dms) == 0:
            raise Exception('No data for selected DM range [%.2f, %.2f],' %
                            (lodm, hidm))

        self.dm2idx = dict([(dm, i) for i, dm in enumerate(self.dms)])
        self.max_downfact = max_downfact

        self.n_success = 0
        self.n_error = 0
        self.n_rejected = 0

        self.dm_delay_map = dict((dm, 0) for dm in self.dms)
        if delays_file:
            # may need more checking (the delays file that is)
            self.dm_delay_map = read_delays(delays_file, self.dm_delay_map)

    def get_t_overlap(self, dm):
        '''
        TBD
        '''
        return self.max_downfact * self.md_map[dm].bin_width + EPSILON

    def is_ok(self, t, dm):
        '''
        Use candidate a (t, dm) yes or no.

        Note: override this in a subclass to implement filtering when reading
        raw candidates.
        '''
        return True

    def iterate_trial(self, dm):
        '''
        Iterate over the candidates in a trial DM with value <dm>.
        '''
        epsilon = EPSILON
        max_n_candidates = HARD_LIMIT_N_CANDIDATES
        sp_file = os.path.join(self.sp_dir, self.sp_map[dm])
        bin_width = self.md_map[dm].bin_width
        i = self.dm2idx[dm]
        dm_delay = self.dm_delay_map[dm]

        with open(sp_file, 'r') as f:
            for line in f:
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
                    if self.is_ok(t, dm):
                        cand = (dm, snr, t, sample, downfact, t - delta_t,
                                t + delta_t, i)

                        self.n_success += 1
                        if self.n_success > max_n_candidates:
                            raise TooManyCandidates(max_n_candidates)
                        yield cand
                    else:
                        self.n_rejected += 1

# --------------------------------------------------------------------------
# -- Write out PRESTO style .singlepulse files -----------------------------


def write_singlepulse_file(filename, candidates):
    candidates.sort(key=lambda x: x[2])
    with open(filename, 'w') as f:
        f.write('''# DM      Sigma      Time (s)     Sample    Downfact\n''')
        for c in candidates:
            f.write('%7.2f %7.2f %13.6f %10d     %3d\n' % (c[0:5]))
