'''
Module contains grouping algorithm that extracts pulses from single pulse data.
'''
from __future__ import division

import copy

import candidate

# --------------------------------------------------------------------------
# -- Pulse class and friends -----------------------------------------------


class Pulse(object):
    def __init__(self, cand, dms_adjacent):
        '''
        Holds a group of single pulse candidate detections that form a pulse.
        '''
        self.head = set([cand])
        self.body = []
        self.interval = [cand[5], cand[6]]
        self.dms_adjacent = dms_adjacent

    def add(self, cand):
        '''
        Add a candidate detection to this pulse.

        Call only when we have a candidate in the pulse for which:
        cand[DM] - c[DM] <= DMS_ADJACENT
        (Where DM is the index of the integer DM entry in the candidate
        tuple.)
        '''
        self.interval = [
            min(self.interval[0], cand[5]),
            max(self.interval[1], cand[6])
        ]
        self.head.add(cand)

    def intersects(self, cand):
        '''
        Returns True when candidate intersects this Pulse.

        Call only with during the construction of pulses with candidates
        from the then highest trial-DM.
        '''

        # Coming in we know that cand[DM] >= c[DM] for
        # any candidate in this candidate pulse. (As that
        # is how this function is called).
        tmp = set()
        ret_val = False
        for c in self.head:
            if cand[7] - c[7] < self.dms_adjacent:
                if cand[5] <= c[6] and c[5] <= cand[6]:
                    ret_val = True
                    break
            else:
                # Never look at c again, only safe if intersects() is called
                # with candidates from the highest trial DM so far (else we
                # will disregard candidates that might intersect other
                # candidates later).
                tmp.add(c)
        self.head = self.head - tmp
        self.body.extend(list(tmp))
        return ret_val

    def merge(self, other):
        '''
        Merge 2 Pulses.

        Call only with Pulses at the same trial-DM.
        '''
        self.interval = [
            min(self.interval[0], other.interval[0]),
            max(self.interval[1], other.interval[1])
        ]

        self.head = self.head | other.head
        self.body.extend(other.body)


def find_peak(pulse):
    '''
    Find highest SNR candidate of pulse.
    '''
    best_candidate = pulse.body[0]
    for c in pulse.body[1:]:
        if c[1] > best_candidate[1]:
            best_candidate = c
    return best_candidate


def find_dm_range(pulse):
    '''
    Determine the range of trial-DMs covered by a pulse.
    '''
    # TODO : Merge with find_peak so that we don't iterate twice
    # (also think about doing this in grouping algo)
    lo_c = pulse.body[0]
    hi_c = pulse.body[0]
    for c in pulse.body[1:]:
        if c[-1] < lo_c[-1]:
            lo_c = c
        elif c[-1] > hi_c[-1]:
            hi_c = c
    return lo_c[0], hi_c[0]


def get_dms(pulse):
    dms = set()

    for c in pulse.body:
        dms.add(c[0])

    return dms


def annotate_pulses(spr, pulses):
    '''
    For each single pulse determine peak snr, dm and dm index.
    '''
    for pulse in pulses:
        best_c = find_peak(pulse)

        pulse.dm = best_c[0]
        pulse.snr = best_c[1]
        pulse.best_i = best_c[7]

        pulse.lo_dm, pulse.hi_dm = find_dm_range(pulse)

    return pulses

# --------------------------------------------------------------------------
# -- Grouping algorithm implementation -------------------------------------


def remove_done(pulses, done):
    '''
    Remove Pulses that can no longer match because of DM value.
    '''
    remove = []
    # For now; clean up pulses that won't grow anymore because their
    # maximum DMs are too low (as evidenced by the fact that there
    # are no more candidates in its ``head'').
    for ii, pulse in enumerate(pulses):
        if not pulse.head:
            remove.append(ii)
    for ii in reversed(remove):
        done.append(pulses.pop(ii))

    return pulses, done


def group(spr, dms_adjacent):
    '''
    Drive the grouping algorithm by calling it for each consecutive trial DM.

    Args:
        spr: A SinglePulseReader (subclass) instance.
        dms_adjacent: An integer; the number of trial DMs below the current
            DM to consider while extracting pulses from candidates.

    Returns:
        A list of extracted pulses.
    '''

    # Note: each Pulse instance has a 'head' against which matches can be
    # made and a 'body' that contains the rest of the candidates. When a
    # candidate is ready its head is empty and all its constituent candidates
    # are in its body. During the candidate matching candidates are removed
    # from the head if their DM is too low to still match (this is safe
    # because the grouping algorithm only moves up in DM). This trick keeps
    # the number of intersections down and makes it easy to recognize pulses
    # that are ready (since their heads are empty).

    pulses = []  # Pulses in the process of being extracted.
    done = []    # List of 'finished' i.e. extracted pulses.

    # Note: spr.dms is a sorted (small->big) list of all trial DMs.
    for dm in spr.dms:
        # Find the appropriate amount of overlap in time for this trial DM.
        t_overlap = spr.get_t_overlap(dm)
        # Add this DM's candidates to the list of pulses. A candidate either
        # gets added to an existing pulse if a match is found or added as a
        # new pulse is no match was found.
        pulses = group_dm(pulses, spr.iterate_trial(dm), t_overlap,
                          dms_adjacent)
        # Remove pulses from the list that can no longer match (because their
        # highest DM candidate is too low to still match any of the new
        # candidates, i.e. whose heads are empty).
        pulses, done = remove_done(pulses, done)

    # For the pulses that were not yet 'done' now that we have reached the
    # highest DM, add them to the list of pulses that are done. (Making
    # sure to empty their heads and adding the candidates to their bodies).
    for p in pulses:
        p.body.extend(list(p.head))
        p.head = set()
    done.extend(pulses)

    print 'There are %d pulses' % len(done)
    print 'Read %d candidates.' % spr.n_success
    return done


def group_dm(pulses, it, t_overlap, dms_adjacent):
    '''
    Group this trial-DM's candidates.

    Note assumes the candidates are sorted in time. XXX (check, t or t - width)
    '''
    min_idx = 0
    tmp_idx = 0

    for cand in it:
        MAXDT2 = 2 * t_overlap
        insert_idx = 0
        matches = []
        for i, pulse in enumerate(pulses[min_idx:]):
            idx = i + min_idx

            # Find valid insertion point:
            if cand[5] > pulse.interval[0]:
                insert_idx = idx + 1

            # Avoid expensive intersection test and find an appropriate
            # place to start the search window for next ``cand'':
            if cand[5] - MAXDT2 > pulse.interval[1]:
                tmp_idx = idx
            elif pulse.intersects(cand):
                matches.append(idx)

            # Cut off the search, no more matches can be found.
            if pulse.interval[0] - MAXDT2 > cand[6]:
                break

        if not matches:
            # No existing pulses match cand, add a new pulse with cand as
            # its only member:
            pulses.insert(insert_idx, Pulse(cand, dms_adjacent))
        else:
            # Matching pulses(s) found, merge them, then add the new cand
            # and finally check the position in the list.
            for ii in reversed(matches[1:]):
                pulses[matches[0]].merge(pulses.pop(ii))
            pulses[matches[0]].add(cand)
            if matches[0] > insert_idx:
                tmp = pulses.pop(matches[0])
                pulses.insert(insert_idx, tmp)

        min_idx = tmp_idx

    return pulses

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


def write_arrivals(filename, pulses):
    '''
    Write out the peak candidates for each of the pulses in a detecion.
    '''

    pulse_peaks = []
    for pulse in pulses:
        pulse_peaks.append(copy.copy(find_peak(pulse)))

    pulse_peaks.sort(key=lambda c: c[2])
    candidate.write_singlepulse_file(filename, pulse_peaks)
