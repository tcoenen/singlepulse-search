'''
Module with pulse train extracting functionality.
'''


def split_bright_dim_pulses(pulses, min_snr, min_n_candidates, lo_dm, hi_dm):
    '''
    Split extracted pulses in above and below threshold.

    Note: consumes the list of Pulses ``pulses''
    '''
    # TODO : rename and refactor this ...

    bright_pulses = []
    dim_pulses = []

    # Throw out anything peaking below DM threshold.
    tmp = []
    if lo_dm is not None:
        while pulses:
            p = pulses.pop(0)
            if p.dm < lo_dm:
                dim_pulses.append(p)
            else:
                tmp.append(p)
        pulses = tmp

    # and peaking above DM threshold:
    tmp = []
    if hi_dm is not None:
        while pulses:
            p = pulses.pop(0)
            if p.dm > hi_dm:
                dim_pulses.append(p)
            else:
                tmp.append(p)
        pulses = tmp

    # do the actual splits
    while pulses:
        p = pulses.pop()
        if p.snr >= min_snr and len(p.body) >= min_n_candidates:
            bright_pulses.append(p)
        else:
            dim_pulses.append(p)

    return bright_pulses, dim_pulses


def extract_pulsetrains(spr, pulses, min_snr, min_n_repeats=0,
                        min_n_candidates=7, n_dms_bright=6, lo_dm=None,
                        hi_dm=None):
    '''
    Find bright pulses that arrive at the same DM.

    Return a list of lists of bright pulses (=Pulse instances)
    and a list of rejects (=Pulse).
    '''
    print 'Looking for DMs with bright single pulses.'
    keepers, rejects = split_bright_dim_pulses(pulses, min_snr,
                                               min_n_candidates, lo_dm, hi_dm)
    keepers.sort(key=lambda x: -x.snr)

    out = []
    while keepers:
        current = keepers.pop(0)

        min_dm = spr.dms[max(0, current.best_i - n_dms_bright)]
        max_dm = spr.dms[min(len(spr.dms) - 1, current.best_i + n_dms_bright)]

        tmp = [current]
        # See whether there are any single pulses above the threshold match the
        # peak DM of the `current' single pulse, remove them from the list of
        # keepers.
        matches = []
        for i, p in enumerate(keepers):
            if min_dm <= p.dm <= max_dm:
                matches.append(i)
        for i in reversed(matches):
            tmp.append(keepers.pop(i))
        if len(tmp) < min_n_repeats:
            # Not enough repeated bright pulse detections for this DM.
            rejects.extend(tmp)
            continue
        out.append(tmp)

    return out, rejects


def add_dim_pulses(spr, bright_pulses, dim_pulses, min_snr=5.5,
                   min_n_candidates=5, n_dms_bright=6):
    '''
    Add low significance pulses to detected bright ones.
    '''
    # bright_pulses is a list of lists of Pulse instances!
    out = []
    while bright_pulses:
        tmp = bright_pulses.pop(0)  # tmp is a list of Pulse instances

        min_dm = spr.dms[max(0, tmp[0].best_i - n_dms_bright)]
        max_dm = spr.dms[min(len(spr.dms) - 1, tmp[0].best_i + n_dms_bright)]
        # See whether there are any other single pulses below the threshold
        # match the peak DM of the list of bright detections under
        # consideration (that list being ``tmp'').
        matches = []
        for i, p in enumerate(dim_pulses):
            if min_dm <= p.dm <= max_dm and p.snr >= min_snr:
                if len(p.body) < min_n_candidates:
                    continue
                matches.append(i)
        for i in reversed(matches):
            tmp.append(dim_pulses.pop(i))
        # Add the list of matching single pulse detections to the output.
        out.append(tmp)

    return out, dim_pulses
