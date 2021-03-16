import os
import logging
from runhelper import make_inputs_from_file
from runhelper import JobFileNameComposer
import cdump2netcdf
from volcMER import HT2unit
import hysplit
from ensemble_tools import ATLtimeloop

logger = logging.getLogger(__name__)

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

def gefs_suffix_list():
    base0 = 'gec'
    base = 'gep'
    suffix = []
    for num in range(0,31,1):
        if num==0:
           suffix.append(base0 + str(num).zfill(2))
        else:
           suffix.append(base + str(num).zfill(2))
    return suffix

def get_cdump_xra_ens(inp, jobid):
    blist = []

    filelocator = JobFileNameComposer(inp["WORK_DIR"], jobid,
                                      inp["jobname"])

    ens_suffix_list = gefs_suffix_list()
    def make_tuple(inval):
        source_tag = "Line to {:1.0f} km".format(inp["top"] / 1000.0)
        suffix = inval[1]
        iii = inval[0] + 1
        cdumpname = "{}.{:03d}".format(
            filelocator.get_cdump_base(stage=iii), iii
        )
        cdumpname = os.path.join(inp["WORK_DIR"], filelocator.get_cdump_filename(stage=0))
        met_tag = suffix
        logger.info("adding to netcdf file :{} {}".format(met_tag, cdumpname))
        return (cdumpname, source_tag, met_tag)

    blist = [make_tuple(x) for x in enumerate(ens_suffix_list)]
    century = 100 * (int(inp["start_date"].year / 100))
    cdumpxra = hysplit.combine_dataset(blist, century=century)
    if cdumpxra.size <= 1:
        logger.debug("ENSEMBLE xra is empty")
    else:
        logger.debug("ENSEMBLE xra is full")
    return cdumpxra

def get_cdump_xra(inp,jobid):
    filelocator = JobFileNameComposer(inp["WORK_DIR"], jobid,
                                      inp["jobname"])
    blist = []
    cdumpname = os.path.join(inp["WORK_DIR"], filelocator.get_cdump_filename(stage=0))
    
    source_tag = "Line to {:1.0f} km".format(inp["top"] / 1000.0)
    met_tag = inp["meteorologicalData"]
    blist = [(cdumpname, source_tag, met_tag)]
    century = 100 * (int(inp["start_date"].year / 100))
    cdumpxra = hysplit.combine_dataset(blist, century=century)
    return cdumpxra

def get_cxra(input_file,wdir,jobid):
    inp = make_inputs_from_file(wdir, input_file).inp
    fname = "xrfile.{}.nc".format(jobid)
    mult = get_conc_multiplier(inp)
    if os.path.isfile(fname):
        logger.info("netcdf file exists. Opening {}".format(fname))
        cxra = xr.open_dataset(fname)
        cxra = cxra.__xarray_dataarray_variable__

        pmult = cxra.attrs["mult"]
        cxra = mult * cxra / pmult
        # remove netcdf file so new one can be written.
        Helper.remove(os.path.join(inp["WORK_DIR"], fname))
    else:
        logger.info("netcdf file does not exist. Creating {}".format(fname))
        cxra = mult * get_cdump_xra_ens(inp,jobid)
    cxra = cxra.assign_attrs({"mult": mult})
    logger.info("writing nc file {}".format(fname))
    cxra = cxra.assign_attrs(inp2attr(inp))
    return cxra

def make_awips_netcdf(input_file,wdir='./',jobid=1002):
    # convert to mg /m3
    cxra = get_cxra(input_file,wdir,jobid)
    inp = make_inputs_from_file(input_file).inp

    # if empty then return an emtpy list.
    if not cxra.size > 1:
        logger.info("make_awips_netcdf: cxra empty. cannot create awips\
                     files")
        return []
    ghash = {}
    ghash["source_latitude"] = inp["latitude"]
    ghash["source_longitude"] = inp["longitude"]
    ghash["source_name"] = inp["VolcanoName"]
    ghash["emission_start"] = inp["start_date"]
    ghash["emission_duration_hours"] = inp["emissionHours"]
    # in mg
    mer = mult / 1e6 / 3600  # kg released in a second.
    ghash["MER"] = mer
    ghash["MER_unit"] = "kg/s"

    #TO DO.
    filenamelocator = JobFileNameComposer(inp["WORK_DIR"],jobid,inp["jobname"])
    awipsname = filenamelocator.get_awips_filename(stage=0)

    c2n = cdump2netcdf.Cdump2Awips(
        cxra, awipsname, munit="mg", jobid=jobid, globalhash=ghash
    )

    awips_files = c2n.create_all_files()
    # returns list of awips files that were created.

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


def example_relative_frequency(input_file='ash_config.txt',wdir='./',jobid=1002):
    # convert to mg /m3
    cxra = get_cxra(input_file,wdir=wdir,jobid=jobid)
    wdir='./'
    inp = make_inputs_from_file(wdir=wdir, config_file = input_file).inp

    # creates ensemble relative frequency plots.
    outfilename  = 'filename.zzz.png'

    # latitude and longitude of volcano.
    vlist = (inp["longitude"],inp["latitude"])

    zlevels = cxra.z.values
    enslist = cxra.ens.values

    # probability levels for plotting.
    clevels = [5,20,40,60,80,95]

    # threshold of exceedance in mg/m3.
    thresh = 0.2

    title = "HYSPLIT ensemble relative frequency exceeding (:0.2f)mg/m3".format(thresh) 
    title += "\n GEFS {} members".format(len(enslist))

    # adjust make sure that ATL plots are not all empty.
    # if maximum value below threshold then adjust threshold so
    # it is 1/10th the max value. some time periods may still be empty.
    adjust = 10

    fignamelist = ATLtimeloop(cxra, enslist, thresh, zlevels,vlist,
                              name = outfilename,
                              norm = True,
                              clevels = clevels,
                              title = title,
                              adjust = adjust)


