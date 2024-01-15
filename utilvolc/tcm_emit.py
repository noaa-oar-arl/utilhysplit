import datetime
import logging
import os
from utilhysplit import emitimes
import pandas as pd

logger = logging.getLogger(__name__)


def construct_efile(
    vals,
    vlat,
    vlon,
    area=1,
    emis_threshold=50,
    vres=1000,
    name="emit.txt",
    date_cutoff=None,
    phash={1:'PASH'},
):
    """
    vals : output from make_outdat method in Inverse class.
    vlat : latitude of vent
    vlon : longitude of vent
    area : area to place ash over.
    emis_threshold : do not use emissions below this value.
    vres : vertical resolution. Needed to create last point in vertical line source.
    phash : dictionary mapping the value in vals (such as 'p006', 'p200') to
            the particle number which is an integer (1,2,3).

    """

    # print("construct_efile function")
    #efile = emitimes.EmiTimes(species=phash)
    efile = emitimes.EmiTimes(species=phash)
    # print("splist", efile.splist)
    # efile.set_species(phash)
    duration = "0100"
    cycle_list = []

    addso4 = False
    if "pSO4" in phash.keys():
        addso4 = True

    for iii, value in enumerate(vals):
        time = pd.to_datetime(value[0])
        if date_cutoff:
            if time >= date_cutoff:
                break
        height = value[1]
        emis = value[2]
        if len(value) > 3:
            part = phash[int(value[3])]
        else:
            part = 1
        # print("construct_efile function:", value, part)
        if emis < emis_threshold:
            emis = 0
        newcycle = False
        if time not in cycle_list:
            newcycle = True
        if newcycle:
            print('adding cycle {}'.format(time))
            efile.add_cycle(pd.to_datetime(time), duration=1)
            cycle_list.append(time)
        print('adding record {:2.2e} {} {}'.format(emis, height, time))
        efile.add_record(
            time,
            duration=duration,
            lat=vlat,
            lon=vlon,
            height=height,
            rate=emis,
            area=area,
            heat=0,
            spnum=part,
            nanvalue=0,
        )
        # SO4 starts with 0 emissions.
        if addso4:
            efile.add_record(
                time,
                duration=duration,
                lat=vlat,
                lon=vlon,
                height=height,
                rate=0,
                area=0,
                heat=0,
                spnum=part,
                nanvalue=0,
            )
        # if there will be a new cycle
        # need to add another record for top height of last point.
        # TO DO. do this for each particle number
        if iii + 1 < len(vals):
            next_time = vals[iii + 1][0]
        else:
            next_time = vals[0][0]
        if next_time != time and vres > 0:
            # do this for each particle number
            efile.add_record(
                time,
                duration=duration,
                lat=vlat,
                lon=vlon,
                height=height + vres,
                rate=0,
                area=area,
                heat=0,
                spnum=part,
                nanvalue=0,
            )
            if addso4:
                efile.add_record(
                    time,
                    duration=duration,
                    lat=vlat,
                    lon=vlon,
                    height=height + vres,
                    rate=0,
                    area=0,
                    heat=0,
                    spnum=part,
                    nanvalue=0,
                )
    # print('writing efile {}', name)
    # efile.write_new(name)
    return efile


def make_efile(
    vals,
    vlat,
    vlon,
    area=1,
    emis_threshold=50,
    vres=1000,
    name="emit.txt",
    date_cutoff=None,
    phash={1:'P060'},
):
    efile = construct_efile(
        vals, vlat, vlon, area, emis_threshold, vres, name, date_cutoff, phash
    )
    efile.write_new(name)
    return efile

