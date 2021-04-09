import os
import logging
import hysplit
import ensemble_tools
from cdump2netcdfAsh import Cdump2Awips
from volcMER import HT2unit

logger = logging.getLogger(__name__)

#----------------------------------------------------------------------------    

def example():
    inp = {}
    # name of output files for ensemble relative frequency plots.
    # zzz will be replaced by plot number for multiple time periods.
    flin = 'douglas_ensrelfreq_zzz.png'
    # name of output files for ensemble mean mass loading plots.
    # zzz will be replaced by plot number for multiple time periods.
    meanflin = 'douglas_ensmeanmass_zzz.png'
    # name for awips2 files
    awipsname = 'out.nc'
 
    inp['VolcanoName'] = 'Douglas' # name of volcano 
    inp['meteorologicalData'] = 'GEFS' # met data set.
    inp['emissionHours'] = 12 # met data set.
    inp['start_date'] = datetime.datetime(2021,4,7,21,15)

    inp['latitude'] = 58.855   # latitude of vent.
    inp['longitude'] = -153.54 # longitude of vent
    inp['eflag'] = 0           # increase(negative) or decrease(positive) MER
    inp['bottom'] = 7021       # altitude of vent (m)
    inp['top'] = 20000         # altitude of plume height (m).
    tdir = './'                # directory to find cdump files in.
    # list of cdump files to use
    cdumplist = ["ashtest_1001_cdump.001",
                 "ashtest_1001_cdump.002",
                 "ashtest_1001_cdump.003",
                 "ashtest_1001_cdump.004"]
    # list of ensemble names corresponding to each cdump file.
    enslist = ["gec00","gep01","gep02","gep03"]
    if len(enslist) != len(cdumplist): 
       print('WARNING: enslist and cdumplist should be same length')
   
    # load data from cdump files into xarray. 
    cxra = maketestra(tdir, cdumplist,enslist)

    # convert unit mass/m3 to mg/m3.
    mult = get_conc_multiplier(inp)   
    cxra = mult * cxra
    inp['mult'] = mult
    cxra = cxra.assign_attrs(inp2attr(inp))

    # create png plots
    plot_ash_rel_freq(cxra,inp,flin,meanflin)

    # create awips2 files
    # dictionary information in ghash will be written into global
    # attributes in the netcdf file.
    ghash["source_latitude"] = inp['latitude']
    ghash["source_longitude"] = inp['longitude']
    ghash['source_name'] = inp['VolcanoName']
    ghash['emission_start'] = inp['start_date']
    ghash['emission_duration_hours'] = inp['emissionHours']
    ghash['MER'] = mult / 1e6 / 3600.0
    ghash['MER_unit'] = 'kg/s'
    c2n = Cdump2Awips(cxra,awipsname,munit='mg',globalhash=ghash)
    fnames = c2n.create_all_files()

 
def maketestblist(dname,cdumplist,enslist):
    # Need list of tuples. (filename, sourcetag, mettag)
    blist = []
    for fname in zip(cdumplist,enslist):
        blist.append((os.path.join(dname, fname[0]), "S1", fname[1]))
    return blist

def maketestra(wdir,cdumplist,enslist):
    blist = maketestblist(wdir,cdumplist,enslist)
    # xrash is an xarray dataset which can be input into
    # Cdump2Awips class initialization.
    xrash = hysplit.combine_dataset(blist,century=2000)
    return xrash

#----------------------------------------------------------------------------    
#----------------------------------------------------------------------------    
#----------------------------------------------------------------------------    
#----------------------------------------------------------------------------    
def inp2attr(inp):
    """
    return dictionary where all values are strings.
    This is so it can be written as attributes to netcdf file.
    """
    atthash = {}
    for key in inp.keys():
        try:
            val = str(inp[key])
        except:
            val = "skip"
        if val != "skip":
            atthash[key] = val
    return atthash


def get_conc_multiplier(inp):
    """
    factor to convert concentration in unit mass /m3 to mg/m3.
    Model must have been run with 1 unit mass/h released.
    """
    height = (inp["top"] - inp["bottom"]) / 1000.0
    # HT2unit outputs in grams. Then multiply by 1e3 to get mg
    conc_multiplier = 1e3 * HT2unit(height) * get_ash_reduction(inp)
    return conc_multiplier

def get_ash_reduction(inp):
    eflag = float(inp["eflag"])
    M63 = 0.01  # fraction of  fine ash
    conc_multiplier = 10 ** (-1 * eflag)
    return M63 * conc_multiplier

def plot_ash_rel_freq(cxra, inp, flin, meanflin):
    """
    plots probability of exceedances/ ensemble relative frequency of exceedance.
    Also plots ensemble mean.
    cxra : xrarrry data array. concentrations should be in units of mg/m3.
    inp  : dictionary with 
           'latitude' and 'longitude' indicating location of volcano.
           'bottom' and 'top' indicating vent height(bottom) and plume height(top)
           'eflag' ash reduction indicator. Positive numbers reduce MER. Negative numbers increase MER.
    flin : str : filename for relative frequency of exceedance plots.
    meanflin : str : filename for ensemble mean plots.
    """
    # Make sure that ATL plots are not all empty.
    # if maximum value below threshold then adjust
    # threshold so it is 1/10 the max value.
    # some time periods may be empty as the adjustment is
    # applied to the entire array.
    adjust = 10

    # must be completed after self.cxra is filled.
    if cxra.size <= 1:
        logger.info("plot_ATL cxra is empty")
        return False
    enslist = cxra.ens.values
    level = cxra.z.values
    thresh = 0.2
    # location of volcano
    vlist = [inp["longitude"], inp["latitude"]]
    logger.debug("NEW FILENAME{}".format(flin))
    clevels = [5, 20, 40, 60, 80, 95]
    title = "HYSPLIT ensemble relative frequency exceeding thresh mg/m3"
    title += "\n GEFS {} members".format(len(enslist))
    fignamelist = ensemble_tools.ATLtimeloop(
        cxra,
        enslist,
        thresh,
        level,
        vlist,
        name=flin,
        norm=True,
        clevels=clevels,
        title=title,
        adjust=adjust,
    )
    # converstion to g/m2 is done in the massload_plot.
    # cxra should still be in units of mg/m3.
    ensemble_tools.massload_plot(cxra, 
                                 enslist,
                                 name=meanflin,
                                 vlist=vlist)
