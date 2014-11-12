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

from ssps import inf

from ssps.candidate.base import SinglePulseReaderMixin
from ssps.candidate.error import DelaysFileWrong, TooManyCandidates
from ssps.candidate.constants import HARD_LIMIT_N_CANDIDATES, EPSILON

# --------------------------------------------------------------------------
# -- Read raw candidates ---------------------------------------------------

SP_PATTERN = re.compile(r'\S+_DM(?P<dm>\d+\.\d+)\.singlepulse')
INF_PATTERN = re.compile(r'\S+_DM(?P<dm>\d+\.\d+)\.inf')

class VanillaBaseReader(object):
    def __init__(self, directory, tstart, tend, delays_file, lodm=None,
                 hidm=None, max_downfact=30):
        '''
        Read PRESTO single pulse data set.

        Note: current implementation assumes the .inf files are in a
        subdirectory of <directory> called INF and the .singlepulse files
        in a subdirectory of <directory> callsed SINGELPULSE.
        '''
        # Vanilla PRESTO single pulse directory handling
        self.sp_dir = os.path.join(directory, 'SINGLEPULSE')
        self.inf_dir = os.path.join(directory, 'INF')
        self.sp_map = self.find_files(self.sp_dir, SP_PATTERN)
        self.inf_map = self.find_files(self.inf_dir, INF_PATTERN)

        # Store desired DM, time range.
        self.tstart = tstart
        self.tend = tend
        self.max_downfact = max_downfact

        # determine trial-DMs, load meta data and extra delays per DM
        dms = list(set(self.sp_map.keys()) & set(self.inf_map.keys()))
        if len(dms) == 0:
            raise Exception('No matching .inf AND .singlepulse files in %s' %
                            directory)

        self.dms = self.grab_dms(dms, lodm, hidm)
        self.md_map = self.read_metadata(self.inf_dir, self.inf_map)
        self.dm_delay_map = self.grab_dm_delay_map(self.dms, delays_file)

        self.dm2idx = dict([(dm, i) for i, dm in enumerate(self.dms)])

        # for book-keeping purposes
        self.n_success = 0
        self.n_error = 0
        self.n_rejected = 0

    def find_files(self, directory, pattern):
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


    def read_metadata(self, directory, dm2inf):
        '''
        Read all given PRESTO style .inf files in directory.
        '''

        metadata_map = {}
        for dm, f in dm2inf.iteritems():
            metadata_map[dm] = inf.inf_reader(os.path.join(directory, f))

        return metadata_map


class SinglePulseReaderCondensed(VanillaBaseReader, SinglePulseReaderMixin):
    def iterate_trial(self, dm):
        # for each candidate return: t
        sp_file = os.path.join(self.sp_dir, self.sp_map[dm])

        with open(sp_file, 'r') as f:
            delay = self.dm_delay_map[dm]
            for line in f:
                split_line = line.split()
                try:
                    t = float(split_line[2])
                    snr = float(split_line[1])
                except:
                    self.n_error += 1
                else:
                    self.n_success += 1
                    if self.tstart <= t + delay <= self.tend:
                        yield t + delay, snr


class SinglePulseReaderBase(SinglePulseReaderMixin, VanillaBaseReader):
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

