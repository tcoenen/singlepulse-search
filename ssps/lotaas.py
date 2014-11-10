# LOTAAS specific format handling
from ssps import inf
from ssps import candidate
import copy

HARD_LIMIT_N_CANDIDATES = 1000000  # Maximum number of candidates in data.
EPSILON = 1e-5  # (Smaller than DM0.00 time-series bin size in seconds!)

class LOTAASBaseReader(object):
    def __init__(self, basename, tstart, tend, delays_file, lodm=None, hidm=None,
                 max_downfact=30):
        '''
        LOTAAS merged single pulse reader metadata handling
        '''
        # construct relevant filenames
        self.singlepulse_file = basename + '.singlepulse'
        self.inf_file = basename + '.inf'

        # scan the data for chunks belonging to 1 DM
        index, binwidth_map = self.scan_data()
        self.index = index

        self.tstart = tstart
        self.tend = tend
        self.max_downfact = max_downfact

        # determine trial-DMs, load meta data and extra delays per DM
        self.dms = self.grab_dms(index.keys(), lodm, hidm)
        self.md_map = self.grab_metadata_map(self.dms, self.inf_file, binwidth_map)
        self.dm_delay_map = self.grab_dm_delay_map(self.dms, delays_file)

        # for book-keeping purposes
        self.n_success = 0
        self.n_error = 0
        self.n_rejected = 0

    def grab_dm_delay_map(self, dms, delays_file):
        '''
        Extract the delays per DM if available.
        '''
        dm_delay_map = dict((dm, 0) for dm in dms)
        if delays_file:
            print 'Loading delays from %s' % delays_file
            dm_delay_map = candidate.read_delays(delays_file, dm_delay_map)

        return dm_delay_map

    def grab_metadata_map(self, dms, inf_file, binwidth_map):
        '''
        Grab the metadata for all intersting DMs.
        '''
        metadata = inf.inf_reader(inf_file)
        metadata_map = {}
        for dm in dms:
            metadata_map[dm] = copy.deepcopy(metadata)
            # LOTAAS specific HACK to get around missing .inf files while
            # still having access to a delay for each DM.
            metadata_map[dm].binwidth = binwidth_map[dm] 

        return metadata_map

    def grab_dms(self, dms, lodm, hidm):
        '''
        Grab the list of intersting DMs (in desired DM range).
        '''
        dms.sort()

        if lodm:
            dms = [dm for dm in dms if lodm <= dm]
        if hidm:
            dms = [dm for dm in dms if dm <= hidm]

        return dms

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


class LOTAASGrabberMixin(LOTAASBaseReader):
    def __init__(self, *args, **kwargs):
        # blah ..
        super(LOTAASGrabberMixin, self).__init__(*args, **kwargs)

        self.dm2idx = dict([(dm, i) for i, dm in enumerate(self.dms)])

    def iterate_trial(self, dm):

        epsilon = EPSILON
        max_n_candidates = HARD_LIMIT_N_CANDIDATES
        bin_width = self.md_map[dm].bin_width
        i = self.dm2idx[dm]
        dm_delay = self.dm_delay_map[dm]

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

    def get_t_overlap(self, dm):
        '''
        tbd
        '''
        return self.max_downfact * self.md_map[dm].bin_width + EPSILON

class LOTAASCondenserMixin(LOTAASBaseReader):
    def iterate_trial(self, dm):

        delay = self.dm_delay_map[dm]
        startpos, endpos = self.index[dm]

        with open(self.singlepulse_file, 'r') as f:
            f.seek(startpos)
            
            while f.tell() < endpos:
                line = f.readline()

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

