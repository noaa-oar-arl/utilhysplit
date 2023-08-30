import datetime
import logging
import os
from subprocess import call

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import seaborn as sns
import xarray as xr

import utilhysplit.hysplit_gridutil as hgu
from monetio.models import hysplit
from utilhysplit import hcontrol
from utilhysplit.evaluation import ensemble_tools, plume_stat
from utilhysplit.metfiles import MetFileFinder
from utilhysplit.plotutils import colormaker
from utilvolc import volcat
from utilvolc.runhelper import Helper

# def save_emis_gefs(inverse):
#    tlist = 'gec00'
#    inverse.save_emis(eii=0, savename = 'gec00')
#    for nlist in np.arange(1,30):
#        ttt = 'gep{:2d}'.format(nlist)
#        inverse.save_emis(eii=nlist

logger = logging.getLogger(__name__)


def check_type_basic(check):
    if isinstance(check, int):
        return True
    if isinstance(check, float):
        return True
    if isinstance(check, str):
        return True
    return False


# def check_type_iterable(check):


def stack_psizes(xrlist, psizelist):
    """
    xrlist : list of data arrays.
    psizelist : list of particle sizes
    lists need to be the same length.
    Returns:
    combined xarray with new ensemble
    dimension which is a tuple of (ens, psize).

    """
    # xrlist = []
    xnewlist = []
    # get list of all xarrays.
    # for xname in namelist[1:]:
    #    xrlist.append(xr.open_dataset(xname))

    # get xbig which is xarray which has domain that
    # encompasses the rest of the domains.
    xbig, xtemp = xr.align(xrlist[0], xrlist[1])
    for iii, xrr in enumerate(xrlist):
        xbig, xtemp = xr.align(xbig, xrr)

    for iii, xrr in enumerate(xrlist):
        if "source" in xrr.coords:
            xrr = xrr.isel(source=0)
        # align to largest domain
        xbig, xnew = xr.align(xbig, xrr)
        # add the psize dimension
        print(psizelist[iii])
        # xnew = xnew.p060
        print(type(xnew))

        xnew = xnew.assign_coords(psize=psizelist[iii]).expand_dims("psize")
        # stack the peize with the ensemble name.
        xnew = xnew.rename({"ens": "ens1"})
        xnew = xnew.stack(ens=("ens1", "psize"))
        xnewlist.append(xnew)
    # concatenate along the ensemble dimension.
    xens = xr.concat(xnewlist, dim="ens")
    return xens


def update_attrs_for_netcdf(dset, datefmt="%Y%m%d %H:%M"):
    for attr in dset.attrs:
        # if attribute is int, float, str then pass.
        check = dset.attrs[attr]
        if check_type_basic(check):
            continue

        if (
            isinstance(check, list)
            or isinstance(check, tuple)
            or isinstance(check, np.ndarray)
        ):
            newlist = []
            for value in check:
                if check_type_basic(value):
                    continue
                if isinstance(value, datetime.datetime):
                    newlist.append(value.strftime(datefmt))
                else:
                    print("unknown type in list", type(check))
                    newlist.append("unknown")
            dset.attrs.update({attr: newlist})
            continue

        if isinstance(value, datetime.datetime):
            dset.attrs.update({attr: value.strftime(datefmt)})
            continue
        print("unknown type ", attr, type(check))
        dset.attrs.update({attr: "unknown"})
    return dset


def testcase(tdir, vdir):
    fname = "xrfile.invbezy.nc"
    bezyloc = [160.587, 55.978]
    vid = "v300250"
    inverse = InverseAsh(tdir, fname, vdir, vid)
    return inverse


def vincent(tdir, vdir):
    # 0.25 degrees
    fname = "xrfile.vincent2.nc"
    # 0.05 degrees
    # fname = 'xrfile.446.nc'
    inverse = InverseAsh(tdir, fname, vdir, vid=None)
    return inverse


def vincent13(tdir, vdir):
    fname = "xrfile.445.nc"
    inverse = InverseAsh(tdir, fname, vdir, vid=None)
    return inverse


# get cdump output.
# pick time period.
# cdump output for that time period.


def getL1L2_filenames(d1, nfiles=100):
    dt = datetime.timedelta(minutes=10)
    flist = []
    dlist = []
    fmt1 = "%Y%j_%H%M"
    d2 = d1
    for iii in np.arange(0, nfiles):
        fname = "geocatL2.GOES-16.Full_Disk.{}30.some_vars.nc".format(d2.strftime(fmt1))
        flist.append(fname)
        dlist.append(d2)
        d2 = d2 + dt
    return flist, dlist


def get_L1L2(tdir, d1, nfiles=100):
    flist, dlist = getL1L2_filenames(d1, nfiles)
    das = []
    for fname in flist:
        fullfname = os.path.join(tdir, fname)
        if os.path.isfile(fullfname):
            das.append(
                volcat.open_dataset(fullfname, decode_times=True, mask_and_scale=False)
            )
        else:
            print("NOT adding {}".format(fullfname))
    return das


def find_area(vmass):
    r2 = vmass.where(vmass > 0)
    r2 = r2.fillna(0)
    r2 = r2.where(r2 <= 0)
    r2 = r2.fillna(1)
    return int(r2.sum())


def plot_L1L2(das, nnn=100):
    masslist = []
    volume = []
    dlist = []
    htlist = []
    alist = []
    for iii in np.arange(0, nnn):
        try:
            vmass = volcat.get_mass(das[iii], clip=True)
        except:
            continue
        vht = volcat.get_height(das[iii], clip=True)
        # print(vmass.sum())
        area = (
            0.03 * 0.03 * 111e3 * 111e3 * np.cos(56 * np.pi / 180.0)
        )  # area of one pixel in meters.
        mass = vmass.sum() * area
        # print('{} kg'.format(mass/1e3))
        # print('{}------'.format(iii))
        masslist.append(float(mass.values))
        dlist.append(das[iii].time.values)
        htlist.append(float(np.max(vht)))
        alist.append(find_area(vmass))
    # print(masslist)
    fig = plt.figure()
    plt.plot(dlist, np.array(masslist) / 1000.0, "--bo")
    ax = plt.gca()
    fig.autofmt_xdate()
    ax.set_ylabel("Total mass (kg)")
    ax.set_xlabel("Time")
    # plt.savefig('bezymass.png')
    plt.show()
    fig = plt.figure()
    plt.plot(dlist, htlist, "--bo")
    ax = plt.gca()
    fig.autofmt_xdate()
    ax.set_ylabel("Maximum height (km)")
    ax.set_xlabel("Time")
    # plt.savefig('bezyheight.png')
    plt.show()
    fig = plt.figue()
    plt.plot(dlist, alist, "--bo")
    ax = plt.gca()
    fig.autofmt_xdate()
    ax.set_ylabel("Number of pixels")
    # plt.savefig('bezypixels.png')
    ax.set_xlabel("Time")


class ParametersIn:
    """
    simple way of changing the n_ctrl, nx_ctrl and
    lbfgs_nbd values in Parmeters_in.dat file.
    """

    def __init__(self, fname):
        # read existing file.
        with open(fname, "r") as fid:
            self.lines = fid.readlines()

    def change_and_write(self, fname, nx_ctrl):
        self.change(nx_ctrl)
        self.write(fname)

    def change(self, nx_ctrl):
        # change relevant lines.
        nx_ctrl = int(nx_ctrl)
        self.lines[1] = " N_ctrl={}\n".format(nx_ctrl)
        self.lines[2] = " Nx_ctrl={}\n".format(nx_ctrl)
        self.lines[40] = " lbfgs_nbd={}*0/\n".format(nx_ctrl)

    def write(self, fname):
        with open(fname, "w") as fid:
            [fid.write(line) for line in self.lines]


def read_emis_df(savename):
    """
    reads dataframe saved by the make_outdat_df method.
    index is height
    columns are date
    values are emissions.
    """
    try:
        emisdf = pd.read_csv(savename, index_col=0, header=1)
    except:
        print("trying this", savename)
        # emisdf = pd.read_csv(savename,index_col='date',header='ht',dtype={'ht':'float'})
        emisdf = pd.read_csv(savename, index_col="date", header="ht")
    return emisdf.dropna()


class InverseOutDat:
    def __init__(self, wdir, fname="out.dat", fname2="out2.dat"):
        # out.dat has estimated release rates in same order as tcm columns.
        # out2.dat has observed(2nd column) and modeled(3rd column) mass loadings.
        self.wdir = wdir
        self.subdir = wdir  # initially set the sub directory to the working directory.

        self.df = pd.DataFrame()  # dataframe which holds data from out.dat
        self.df2 = pd.DataFrame()  # dataframe which holds data from out2.dat

        self.fname = fname  # name of out.dat file
        self.fname2 = fname2  # name of out2.dat file

    def read_out(self, name):
        """
        name : str : filename of out.dat file which has estimated release rates in same
                     order as tcm columns. Currently first column is just a dummy.

        Returns :
        df : pandas dataframe/
        """

        # out.dat has estimated release rates in same order as tcm columns.
        wdir = self.wdir
        if os.path.isfile(os.path.join(wdir, name)):
            df = pd.read_csv(os.path.join(wdir, name), sep="\s+", header=None)
        else:
            print("File not found {}".format(os.path.join(wdir, name)))
            df = pd.DataFrame()
        return df

    def get_emis(self, name=None):
        """
        Reads in out.dat file
        """
        if not name:
            name = self.fname
        else:
            self.fname = name
        df = self.read_out(name)
        self.df = df
        return df

    def emis_hist(self, units="g/h"):
        if units == "g/h":
            vals = np.log10(self.df[1].values)
        nmin = int(np.floor(np.min(vals)))
        nmax = int(np.ceil(np.max(vals)))
        nbins = len(np.arange(nmin, nmax, 1))
        plt.hist(vals, bins=nbins)
        ax = plt.gca()
        ax.set_xlabel("Log emission (g/h)")

    def get_conc(self, name=None):
        """
        Reads in out2.dat file
        """
        if not name:
            name = self.fname2
        df = self.read_out(name)
        df.columns = ["index", "observed", "model"]
        self.df2 = df
        return df

    def plot_conc(self, cmap="viridis"):
        sns.set()
        sns.set_style("whitegrid")
        df = self.df2
        # plt.plot(df['observed'],df['model'],'k.',MarkerSize=3)
        cb = plt.hist2d(
            df["observed"],
            df["model"],
            cmap=cmap,
            norm=mpl.colors.LogNorm(),
            bins=[20, 20],
        )
        cbar = plt.colorbar(cb[3])
        cbar.ax.set_ylabel("Number of Points")
        ax = plt.gca()
        ax.set_xlabel("observed")
        ax.set_ylabel("modeled")
        nval = np.max(df["observed"])
        # plot 1:1 line
        plt.plot([0, nval], [0, nval], "--b", linewidth=1)
        return ax

    def get_vals(self):
        return -1


def get_inp_hash(wdir, configfile):
    from utilvolc.runhelper import make_inputs_from_file

    setup = make_inputs_from_file(wdir, configfile)
    setup.add_inverse_params()
    return setup.inp


def get_sourcehash(wdir, configfile):
    from utilvolc.invhelper import inverse_get_suffix_list
    from utilvolc.runhelper import make_inputs_from_file

    setup = make_inputs_from_file(wdir, configfile)
    setup.add_inverse_params()
    sourcehash = inverse_get_suffix_list(setup.inp)
    return sourcehash


# For the inverse modeling work flow goes
# inverse = ai.InverseAsh(...)
# loop through inverse.perpare_one_time(....)
# use compare_plotsA to check coverage
# inverse.make_tcm_mult() to create TCM
# plot_tcm to check.
# inverse.write_tcm
# Run the TCM

# call to InverseOutDat class to read output from inversion algorithm.
# input dataframe from InverseOutDat into
# inverse.plot_outdat()


def create_runtag(tag, tii, remove_cols, remove_rows, remove_sources, remove_ncs=0):
    base = tag
    times = str.join("_", map(str, tii))
    tag2 = ""
    if remove_cols:
        tag2 += "T"
    else:
        tag2 += "F"
    if remove_rows:
        tag2 += "T"
    else:
        tag2 += "F"
    if remove_ncs > 0:
        tag3 = "w{}".format(remove_ncs)
    else:
        tag3 = ""
    if remove_sources:
        tag4 = "_{}".format(str.join("_", list(map(str, remove_sources))))
    else:
        tag4 = ""
    rval = "Run{}_{}_{}{}{}".format(base, times, tag2, tag3, tag4)
    return rval


class InverseAshEns:
    """
    Inverse runs from a meteorological ensemble.
    May also be used for a single run.
    Utilizes a list of InverseAsh objects each of which represents 1 run.
    """

    def __init__(
        self,
        tdirlist,
        fnamelist,
        vdir,
        vid,
        configdir="./",
        configfile=None,
        ensdim="ens",
        verbose=False,
    ):
        self.fnamelist = fnamelist  # list of netcdf files with HYSPLIT output (cdump)
        self.tcm_names = []  # list of names of tcm files written.
        self.n_ctrl_list = []  # list of numbers for Parameter_in.dat
        self.vdir = vdir  # directory where volcat information is found.
        self.phash = {}  # dictionary with particle size information
        self.taglist = self.create_taglist(fnamelist)

        # create self.invlist with add_invlist method.
        self.invlist = []  # list of  InverseAsh objects
        self.add_invlist(tdirlist, vdir, vid, configdir, configfile, ensdim, verbose)

        # directories
        # use set_directory method to change
        default = "./"
        self.datadir = default  # location of netcdf files with HYSPLIT output (cdump)
        self.wdir = default  # location to write outputs, create new directories
        self.subdir = default  # name of subdirectory in working directory
        self.execdir = default  # directory where inversion code executable is located
        self.hysplitdir = default  # directory where HYSPLIT executables are lcated.

        self.subdirlist = []  # list of subdirectories.
        self.tcm_names = []  # list of names of tcm files written.
        self.n_ctrl_list = []  # list of numbers for Parameter_in.dat

    def create_taglist(self, fnamelist):
        try:
            taglist = [self.make_tag(x) for x in fnamelist]
        except:
            logger.warning(
                "InverseAshEns string not in expected form {}".format(fnamelist[0])
            )
            taglist = list(map(str, np.arange(0, len(fnamelist))))
        return taglist

    @staticmethod
    def make_tag(fname):
        if "gep" in fname or "gec" in fname:
            # These are using GEFS files.
            # assume strings of form NAME_gep04.nc
            # assume strings of form NAME_gec00.nc
            return fname[-8:-3]
        else:
            temp = fname.split(".")
            # assume strings of form xrfile.TAG.nc
            if len(temp) == 3:
                return temp[1]
            # assume strings of form TAG.nc
            elif len(temp) == 2:
                return temp[0]

    def close_arrays(self):
        """
        InverseAshEns class
        """
        for inv in self.invlist:
            inv.close_arrays()

    def add_invlist(self, tdirlist, vdir, vid, configdir, configfile, ensdim, verbose):
        """
        InverseAshEns class
        """
        if len(tdirlist) != len(self.fnamelist):
           logger.warning('Length of directory list not the same as filenamelist')
        for hruns in zip(tdirlist, self.fnamelist):
            self.invlist.append(
                InverseAsh(
                    hruns[0],
                    hruns[1],
                    vdir,
                    vid,
                    configdir,
                    configfile,
                    ensdim=ensdim,
                    verbose=verbose,
                )
            )

    def add_config_info(self, configdir, configfile):
        """
        InverseAshEns class
        """
        for inv in self.invlist:
            inv.add_config_info(configdir, configfile)

    def set_directory(self, wdir, execdir, datadir, hysplitdir, subdir=None):
        """
        InverseAshEns class
        """
        self.datadir = datadir
        self.wdir = wdir
        if not subdir:
            self.subdir = wdir
        else:
            self.subdir = os.path.join(wdir, subdir)
        self.execdir = execdir
        self.hysplitdir = hysplitdir

    def print_directories(self):
        """
        InverseAshEns class
        """
        print("Working directory, wdir :{}".format(self.wdir))
        print("execdir :{}".format(self.execdir))
        print("hysplitdir :{}".format(self.hysplitdir))
        print("subdir :{}".format(self.subdir))
        print("vdir :{}".format(self.vdir))
        print("datadir :{}".format(self.datadir))

    def set_subdirectory(self, runtag):
        """
        InverseAshEns class
        """
        # make subdirectories for different TCM options.
        sdir = os.path.join(self.wdir, runtag)
        if not os.path.isdir(sdir):
            callstr = "mkdir -p " + os.path.join(sdir)
            call(callstr, shell=True)
        self.subdir = sdir
        if sdir not in self.subdirlist:
            self.subdirlist.append(sdir)
        return sdir

    def set_concmult(self, mult):
        """
        InverseAshEns class
        """
        for hrun in self.invlist:
            hrun.set_concmult(mult)

    def reset_tcm(self):
        """
        InverseAshEns class
        """
        self.tcm_names = []  # list of names of tcm files written.
        self.n_ctrl_list = []  # list of numbers for Parameter_in.dat

    def get_mean_tcm(self):
        """
        InverseAshEns class
        Not implemented yet.
        """
        for hrun in zip(self.invlist, self.taglist):
            temp = 1

    def write_tcm(self, tcm_name, reset=True, verbose=False):
        #InverseAshEns 
        """
        """
        if reset:
            self.reset_tcm()
        for hrun in zip(self.invlist, self.taglist):
            tname = tcm_name.replace(".txt", "")
            tname = "{}_{}.txt".format(tname, hrun[1])
            hrun[0].write_tcm(tname)
            if verbose:
                print("writing TCM {} {}".format(os.getcwd(), tname))
            if tname not in self.tcm_names:
                self.tcm_names.append(tname)
                self.n_ctrl_list.append(hrun[0].n_ctrl)

    def plot_tcm(self, ensi=None):
        # InverseAshEns 
        """
        """
        if ensi:
            self.invlist[ensi].plot_tcm()
            return True
        for hrun in self.invlist:
            hrun.plot_tcm()
            plt.show()
        return True

    def make_tcm_names(self):
        #InverseAshEns 
        out_name1 = "out.dat"
        out_name2 = "out2.dat"
        name1 = []
        name2 = []
        for tag in self.taglist:
            name1.append("{}_{}".format(tag, out_name1))
            name2.append("{}_{}".format(tag, out_name2))
        return name1, name2

    def create_emit_output(self, outname, source="M4", overwrite=True, attrs=None):
        #InverseAshEns 
        """
        creates netcdf file with outputs from the runs with the emit-times files.
        This can be fed into the AshEval class.
        attrs : dictionary with any additional information to be added to global attributes.
        """
        print("making file", outname)
        blist = []
        if os.path.isfile(outname):
            if not overwrite:
                print("File already exists {}".format(outname))
                return -1
            else:
                print("Overwriting file {}".format(outname))
                Helper.remove(outname)
        newattr = {
            "sources": list(self.invlist[0].tcm_columns),
            "Method": "lbfgsb",
            "units": "g/m3",
        }
        if attrs:
            newatrr.update(attr)
        for tag in self.taglist:
            newattr["MetData"] = "{}".format(tag)
            blist.append(
                (os.path.join(self.subdir, "cdump.{}".format(tag)), source, tag)
            )
        dset = hysplit.combine_dataset(blist)
        dset.attrs.update(newattr)
        # try:
        #    update_attrs_for_netcdf(dset)
        # except:
        #    print('attr update did not work')
        #    pass
        try:
            dset.to_netcdf(outname)
        except:
            print("netcdf not created")
            pass
        return dset

    def run_hysplit(self):
        #InverseAshEns 
        """
        """
        import time
        from utilhysplit.runhandler import ProcessList

        processhandler = ProcessList()
        processhandler.pipe_stdout()
        processhandler.pipe_stderr()
        os.chdir(self.subdir)
        cmd = os.path.join(self.hysplitdir, "exec", "hycs_std")
        for tag in self.taglist:
            print("running ", cmd, tag, os.getcwd())
            processhandler.startnew([cmd, tag], self.subdir, descrip=tag)
            time.sleep(5)
        #    Helper.execute([cmd,tag])
        done = False
        seconds_to_wait = 30
        total_time = 0
        max_time = 60 * 6000
        while not done:
            num_procs = processhandler.checkprocs()
            print("in loop {}s procs {}".format(total_time, num_procs))
            if num_procs == 0:
                done = True
            time.sleep(seconds_to_wait)
            total_time += seconds_to_wait
            if total_time > max_time:
                processhandler.checkprocs()
                processhandler.killall()
                print("WARNING: HYSPLIT run timed out")
                done = True

    def run_tcm(self):
        #InverseAshEns 
        """
        """
        out_name1 = "out.dat"
        out_name2 = "out2.dat"
        inp_name = "TCM_sum.csv"
        cmd = os.path.join(self.execdir, "new_lbfgsb.x")
        new_name1, new_name2 = self.make_tcm_names()
        for iii, tcm in enumerate(self.tcm_names):
            print("run_tcm tag", self.taglist[iii])
            os.chdir(self.subdir)
            print("working in ", os.getcwd())
            params = ParametersIn(os.path.join(self.wdir, "Parameters_in.dat.original"))
            params.change_and_write(
                os.path.join(self.subdir, "Parameters_in.dat"), self.n_ctrl_list[iii]
            )

            Helper.remove(inp_name)

            # remove any existing output files.
            Helper.remove(out_name1)
            Helper.remove(out_name2)

            # run

            Helper.copy(tcm, inp_name)
            print(cmd)
            # Helper.execute_with_shell(cmd)
            Helper.execute(cmd)
            # move output files to new names
            Helper.move(out_name1, new_name1[iii])
            Helper.move(out_name2, new_name2[iii])
            Helper.move("fort.188", "fort.188.{}".format(self.taglist[iii]))

    def plot_outdat(self, eii=None, name=None):
        # InverseAshEns
        """
        eii : int or list of ints.

        Plot of modeled vs. observed values.

        """
        ilist = self.read_outdat(eii)
        for iii, io in enumerate(ilist):
            if not eii:
                nnn = iii
            else:
                if isinstance(eii, int):
                    nnn = eii
                else:
                    nnn = eii[iii]
            print(self.taglist[nnn])
            df2 = io.get_conc(name=name)
            # df = io.get_emis()
            try:
                ax = io.plot_conc()
            except:
                pass
            plt.show()
            # io.emis_hist()
            # ax = plt.gca()
            # ax.set_xlim(-1,10)
            # plt.show()

    def read_outdat(self, eii=None):
        #InverseAshEns 
        """
        eii : int or list of ints.

        Returns:
        ilist : list of InverseOutDat objects

        """
        # ensemble
        name1, name2 = self.make_tcm_names()
        ilist = []
        if isinstance(eii, int):
            eii = [eii]
        for iii, outdat in enumerate(zip(name1, name2)):
            if not eii or (iii in eii):
                # if eii:
                #    print("read outdat", outdat)
                # print("read outdat", outdat)
                # print(self.subdir, outdat[0])
                io = InverseOutDat(self.subdir, outdat[0], outdat[1])
                ilist.append(io)
        return ilist

    def make_emit_name(self, subdir, tag):
        # InverseAshEns 
        """
        """
        return os.path.join(self.subdir, "{}_{}".format(tag, self.emitname))

    def add_phash(self, phash=None):
        # InverseAshEns 
        if not isinstance(phash, dict):
            self.phash = self.invlist[0].phash
        else:
            self.phash = phash

    def make_efile(self, vloc, emis_threshold=1e5, eii=None):
        # InverseAshEns 
        """
        creates emit-times, CONTROL and SETUP files.
        vloc : [longitude, latitude]
        """
        print("Working directory", os.getcwd())
        self.emitname = "emit.txt"
        ilist = self.read_outdat(eii)
        for iii, io in enumerate(ilist):
            tag = self.taglist[iii]
            emit_name = self.make_emit_name(self.subdir, tag)
            df = io.get_emis()
            if df.empty:
                print("No emissions generated")
                continue
            vals = self.invlist[iii].make_outdat(df)
            area = self.invlist[iii].inp["area"]
            efile = make_efile(
                vals,
                vlat=vloc[1],
                vlon=vloc[0],
                area=area,
                emis_threshold=emis_threshold,
                name=emit_name,
                phash=self.phash,
            )
            inp = self.invlist[0].inp.copy()
            if "gefs" in inp["meteorologicalData"].lower():
                # inp['forecastDirectory'] = os.path.join(self.datadir,'gefsnew')
                print("Forecast directory", inp["forecastDirectory"])
                inp["archivesDirectory"] = os.path.join(self.datadir, "wrong")
                metstr = inp["meteorologicalData"]
            elif "gfs0p25" in inp["meteorologicalData"].lower():
                print("Archive  directory", inp["archivesDirectory"])
                #    inp["archivesDirectory"] = os.path.join(self.datadir, "")
                #    inp["forecastDirectory"] = os.path.join(self.datadir, "")
                metstr = "gfs0p25"
            elif "era5" in inp["meteorologicalData"].lower():
                inp["archivesDirectory"] = os.path.join(self.datadir, "")
                inp["forecastDirectory"] = os.path.join(self.datadir, "")
                metstr = "era5"
            elif "merra2" in inp["meteorologicalData"].lower():
                inp["archivesDirectory"] = os.path.join(self.datadir, "")
                inp["forecastDirectory"] = os.path.join(self.datadir, "")
                metstr = "merra2"
            else:
                inp["archivesDirectory"] = os.path.join(self.datadir, "")
                inp["forecastDirectory"] = os.path.join(self.datadir, "")
                metstr = "unknown"
            make_control(
                efile,
                self.wdir,
                "CONTROL.default",
                self.subdir,
                tag,
                forecast_directory=inp["forecastDirectory"],
                archive_directory=inp["archivesDirectory"],
                lat=float(inp["latitude"]),
                lon=float(inp["longitude"]),
                metstr=metstr,
            )
            make_setup(self.wdir, "SETUP.default", self.subdir, suffix=tag)
            Helper.copy(
                os.path.join(self.wdir, "ASCDATA.CFG"),
                os.path.join(self.subdir, "ASCDATA.CFG"),
            )
            # return efile

    def save_emis(self, savename, eii=None, name=None, ens=True):
        # InverseAshEns 
        """
        Saves emissions in a csv file
        """
        ilist = self.read_outdat(eii)
        basename = savename
        for iii, io in enumerate(ilist):
            sname = io.fname
            sname = sname.replace("_out.dat", "")
            if ens:
                savename = "{}{}".format(sname, basename)
            else:
                savename = basename
            df = io.get_emis(name=name)
            savename = os.path.join(self.subdir, savename)
            self.invlist[0].make_outdat_df(df, savename=savename, part="condense")

    def plot_outdat_all(
        # InverseAshEns 
        self, eii=None, unit="kg/s", profile=False, ax=None, name=None, subdirlist=None
    ):
        """
        For plotting files in more than one subdirectory together.
        """
        if isinstance(subdirlist, str):
            subdirlist = [subdirlist]
        if not isinstance(subdirlist, list):
            subdirlist = self.subdirlist
        cmap = "viridis"
        cm = colormaker.ColorMaker(
            cmap, len(subdirlist), ctype="hex", transparency=None
        )
        cmlist = cm()
        sns.set_style("whitegrid")
        sns.set_context("notebook")
        labels = []
        if not ax:
            fig = plt.figure(1)
            ax = fig.add_subplot(1, 1, 1)
        for iii, sdir in enumerate(subdirlist):
            self.subdir = sdir
            labels.append(sdir.split("/")[-1])
            ax, axg = self.plot_outdat_ts(
                eii=eii,
                unit=unit,
                clr="#" + cmlist[iii],
                profile=profile,
                ax=ax,
                name=name,
            )
        handles, junk = ax.get_legend_handles_labels()
        figlegend = plt.figure(2)
        sns.set_style("white")
        axg = figlegend.add_subplot(1, 1, 1)
        axg.legend(handles, labels, loc="center", fontsize=20)
        axg.axis("off")
        plt.tight_layout()
        return ax, axg
        # plt.savefig('emissions{}_legend.png'.format(tag))
        # plt.show()

    def plot_outdat_ts(
        self, eii=None, unit="kg/s", clr=None, profile=False, ax=None, name=None
    ):
        # InverseAshEns 
        """
        eii :
        """
        # Ensemble.
        dflist = []
        sns.set_style("whitegrid")
        sns.set_context("notebook")
        if not ax:
            fig = plt.figure(1)
            ax = fig.add_subplot(1, 1, 1)
        # list of InverseOutDat objects
        ilist = self.read_outdat(eii)

        if not clr:
            cmap = "viridis"
            cm = colormaker.ColorMaker(cmap, len(ilist), ctype="hex", transparency=None)
            clrlist = cm()
            clrlist = ["#" + x for x in clrlist]
            # clrlist = ["--k", "--r", "--b", "--g", "--c", "--y", "--m"]
        else:
            clrlist = [clr]

        jjj = 0
        for iii, io in enumerate(ilist):
            df = io.get_emis(name=name)
            if profile:
                ax, df = self.invlist[0].plot_outdat_profile(
                    df,
                    ax=ax,
                    clr=clrlist[jjj],
                    unit=unit,
                )
            else:
                ax, df = self.invlist[0].plot_outdat_ts(
                    df, ax=ax, clr=clrlist[jjj], unit=unit
                )
            df2 = df.copy()
            df2.columns = df2.columns.map(lambda x: x[1])
            dflist.append(df2)
            jjj += 1
            if jjj >= len(clrlist):
                jjj = 0
        plt.xticks(rotation=45)

        # for plotting the mean
        temp = pd.concat(dflist)
        temp = temp.reset_index(names='level').groupby('level').mean()
        if profile:
            ax, mass = plot_outdat_profile_function(
                temp, ax=ax, clr="r", label="mean", lw=5, alpha=0.5
            )
        else:
            plot_outdat_ts_function(temp, ax=ax, label='mean', clr="r", lw=5, alpha=0.5)

        handles, labels = ax.get_legend_handles_labels()
        figlegend = plt.figure(2)
        sns.set_style("white")
        axg = figlegend.add_subplot(1, 1, 1)
        axg.legend(handles, labels, loc="center", fontsize=20)
        axg.axis("off")
        plt.tight_layout()

        return ax, axg

    def make_tcm_mult(
        self,
        tiilist,
        remove_cols=True,
        remove_rows=True,
        remove_sources=None,
        remove_ncs=0,
    ):
        # InverseAshEns 
        """
        calls make_tcm_mult for each instance of  InverseAsh in invlist.
        make_tcm_mult creates a TCM for multiple time periods.
        """
        for hrun in self.invlist:
            a, b, c, = hrun.make_tcm_mult(
                tiilist, remove_cols, remove_rows, remove_sources, remove_ncs
            )

    def write_descrip(self, tiilist, remove_cols, remove_rows, remove_sources):
        # InverseAshEns 
        """
        """
        with open(os.path.join(self.subdir, "readme.txt"), "a") as fid:
            fid.write("-----------------")
            fid.write(datetime.datetime.now().strftime("%Y %m %d %HH:%MM UTC \n"))
            fid.write("Dates used in TCM")
            for tii in tiilist:
                date = self.invlist[0].cdump.times.values[tii]
                fid.write("{}  {}".format(tii, date))
                fid.write("Columns Removed")

    def prepare_one_time(self, daterange, st="start", zvals=None, verbose=False):
        # InverseAshEns 
        """
        """

        for iii, hrun in enumerate(self.invlist):
            if hrun.get_sampling_time() != st:
                logger.warning("sampling time does not match")
                hrun.change_sampling_time(st)
            hrun.prepare_one_time(daterange, zvals=zvals)
            try:
                hrun.prepare_one_time(daterange, zvals=zvals, verbose=verbose)
            except Exception as eee:
                logger.warning("ERROR in preparing time {}".format(self.taglist[iii]))
                logger.warning(str(eee))

    def compare_plotsA(self, daterange=None, tii=None, zii=None, vloc=None):
        # InverseAshEns 
        """
        daterange : list of datetime.datetime objects
        tii : int
        """
        clist = []
        for iii, hrun in enumerate(self.invlist):
            print("tag", tii, self.taglist[iii])
            if not daterange:
                daterange = [datetime.datetime.now(), datetime.datetime.now()]
            csum = hrun.compare_plotsA(daterange=daterange, tii=tii, zii=zii)
            clist.append(csum)
            if vloc:
                plt.plot(vloc[0], vloc[1], "y^", MarkerSize=20)
            plt.show()
        # if len(clist)>1:
        #   cxra = xr.concat(clist)
        return clist


class InverseAshPartEns(InverseAshEns):
    """
    Inverse runs with different particle sizes.
    """

    def __init__(
        self,
        tdirlist,
        fnamelist,
        vdir,
        vid,
        configdir="./",
        configfile=None,
        ensdim="ens",
        taglist=["P"],
        verbose=False,
    ):
        self.invlist = []  # list of  InverseAsh objects
        self.fnamelist = fnamelist
        self.tcm_names = []  # list of names of tcm files written.
        self.n_ctrl_list = []  # list of numbers for Parameter_in.dat
        self.vdir = vdir

        self.taglist = taglist

        self.invlist.append(
            InverseAshPart(
                tdirlist,
                fnamelist,
                vdir,
                vid,
                configdir,
                configfile,
                ensdim=ensdim,
                verbose=verbose,
            )
        )

        # for hruns in zip(tdirlist, fnamelist):
        #    self.invlist.append(InverseAsh(hruns[0],hruns[1],vdir,vid,
        #                                   configdir,configfile,
        #                                   ensdim=ensdim,
        #                                   verbose=verbose))


class InverseAsh:
    def __init__(
        self,
        tdir,
        fname,
        vdir,
        vid,
        configdir="./",
        configfile=None,
        verbose=False,
        ensdim="ens",
    ):
        """
        configfile : full path to configuration file.

        fname : str. full file path to netcdf file with model data.
                xarray DataArray or DataSet with cdump information
        """
        # InverseAsh class

        self.tdir = tdir  # directory for hysplit output
        self.fname = fname  # name of hysplit output

        self.vdir = vdir  # directory for volcat data
        self.vid = vid  # volcano id.

        self.n_ctrl = 0  # determined in the write_tcm method.
        # needed for Parameters_in.dat input into inversion algorithm

        # keep volcat arrays for different averaging times.
        self.volcat_hash = {}

        # prepare_one_time method adds data to these dictionaries.
        # the data in cdump_hash and volcat_avg_hash have been aligned.
        # and are on the same grid.
        self.volcat_avg_hash = {}
        self.volcat_ht_hash = {}
        self.cdump_hash = {}

        # multiplication factor if more than 1 unit mass released.
        self.concmult = 1

        # get_cdump method loads the cdump file.
        self.cdump = None
        self.get_cdump(tdir, fname, verbose, ensdim)

        # add_config_info creates these two dictionaries from the configuration.
        self.sourcehash = {}
        self.inp = {}
        self.add_config_info(configdir, configfile)

        # These values are for if spatial coarsening is used.
        self.coarsen = None
        self.original_volcat_avg_hash = {}
        self.original_volcat_ht_hash = {}
        self.original_cdump_hash = {}

        self.set_bias_correction(slope=None, intercept=None)

        # particle size information

        # tcm_columns will be filled in later.
        self.tcm_columns = None

    def close_arrays(self):
        """ """
        # InverseAsh class
        self.cdump.close()
        for key in self.volcat_avg_hash.keys():
            self.volcat_avg_hash[key].close()
        # self.massload.close()

    def copy(self):
        """ """
        # InverseAsh class
        iacopy = InverseAsh(self.tdir, self.fname, self.vdir, self.vid, self.n_ctrl)
        iacopy.volcat_avg_hash = self.volcat_avg_hash
        iacopy.volcat_cdump_hash = self.cdump_hash
        icacopy.concmult = self.concmult
        return icacopy

    def print_summary(self):
        """ """
        # InverseAsh class
        print("Observations availalbe in volcat_avg_hash")
        print(self.volcat_avg_hash.keys())
        print("times in cdump file")
        self.print_times()

    def add_cdump(self,tdir,fname,verbose=False,ensdim='ens'):
        """
        adds another cdump file along the time dimension.
        """
        cdump = self.get_cdump_sub(tdir,fname,verbose,ensdim)
        return xr.concat([self.cdump,cdump],dim='time')

    def get_cdump(self,tdir,fname,verbose=False,ensdim='ens'):
        self.cdump = self.get_cdump_sub(tdir,fname,verbose,ensdim)

    def get_cdump_sub(self, tdir, fname, verbose=False, ensdim="ens"):
        """ """
        # InverseAsh class
        # hysplit output. xarray.
        print("working on {}".format(fname))
        if isinstance(fname, str):
            cdump = xr.open_dataset(os.path.join(tdir, fname), engine="netcdf4")
            logger.info("opening {} {}".format(tdir, fname))
        else:
            cdump = fname
        if not hysplit.check_grid_continuity(cdump):
            logger.info("Grid not continuous")
        if not isinstance(cdump, xr.core.dataarray.DataArray):
            # turn dataset into dataarray
            temp = list(cdump.keys())
            cdump = cdump[temp[0]]
        # get rid of source dimension
        if ensdim == "ens":
            cdump = cdump.isel(source=0)
        elif ensdim == "source":
            cdump = cdump.isel(ens=0)
        cdump, dim = ensemble_tools.preprocess(cdump)
        if dim == "source":
            cdump = cdump.rename({"ens": "metid"})
            cdump = cdump.rename({"source": "ens"})
        return cdump.fillna(0)
        #self.cdump = cdump.fillna(0)
        # print('HERE')
        # self.cdump = cdump

    def add_config_info(self, configdir, configfile):
        """ """
        # InverseAsh class
        # the ens dimension holds is key to what emission source was used.
        # the sourcehash is a dictionary
        # key is the ensemble number
        # values is another dictionary with
        # sdate: begin emission
        # edate: end emission
        # bottom : lower height of emission
        # top : upper height of emission.
        if configfile:
            if not os.path.isfile(os.path.join(configdir, configfile)):
                configfile = None
        if configfile:
            self.sourcehash = get_sourcehash(configdir, configfile)
            self.inp = get_inp_hash(configdir, configfile)
        else:
            self.sourcehash = {}
            self.inp = {}

    def add_inp(self, configdir, configfile):
        """ """
        # InverseAsh class
        self.inp = get_inp_hash(configdir, configfile)

    def set_concmult(self, concmult):
        """ """
        # InverseAsh class
        self.concmult = concmult

    def get_volcat(self, daterange, verbose=False):
        """ """
        # InverseAsh class
        vdir = self.vdir
        vid = self.vid
        tii = self.time_index(daterange[0])
        # if tii not in self.volcat_hash.keys():
        done = True
        if done:
            das = volcat.get_volcat_list(
                vdir,
                daterange=daterange,
                vid=vid,
                decode_times=True,
                verbose=verbose,
                include_last=False,
            )
            self.volcat_hash[tii] = das
        else:
            das = self.volcat_hash[tii]
        # create one dataset with dimension of time.
        if len(das) > 0:
            # vset = xr.concat(das, dim="time")
            vset = volcat.combine_regridded(das)
        else:
            print("No volcat files found ")
            return xr.DataArray()
        return vset

    def clip(self, dummy, buf=0):
        """ """
        # InverseAsh class
        # clip ra according to where dummy has 0's.
        aaa = np.where(dummy > 0)
        a1 = np.min(aaa[0])
        a2 = np.max(aaa[0])
        b1 = np.min(aaa[1])
        b2 = np.max(aaa[1])

        if a2 + buf < dummy.y.values[-1]:
            a2 = a2 + buf

        if b2 + buf < dummy.x.values[-1]:
            b2 = b2 + buf

        a1 = a1 - buf
        b1 = b1 - buf
        if a1 < 0:
            a1 = 0
        if b1 < 0:
            b1 = 0

        return a1, a2, b1, b2

    def print_times(self):
        """ """
        # InverseAsh class
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        for time in timelist:
            print(time.strftime("%Y %m %d %H:%Mz"))

    def get_time(self, tii):
        """
        tii : integer
        """
        # InverseAsh class
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        return timelist[tii]

    def time_index(self, time):
        """ """
        # InverseAsh class
        timelist = [pd.to_datetime(x) for x in self.cdump.time.values]
        try:
            iii = timelist.index(time)
        except:
            # print('timelist', timelist)
            iii = None
        return iii

    def make_tcm_mult(
        self,
        tiilist,
        remove_cols=True,
        remove_rows=True,
        remove_sources=None,
        remove_ncs=0,
    ):
        # InverseAsh class
        """
        Creates a tcm for multiple time periods.

        tiilist : list of integers indicating time periods
        Return:
        t3  : numpy array with tcm matrix
        lat : numpy array with latitudes
        lon : numpy array with longitudes

        """
        # make the tcm for multiple time periods.
        tcmlist = []
        latlist = []
        lonlist = []
        for tii in tiilist:
            print(self.cdump.time.values[tii])
            tcm, model_lat, model_lon, columns = self.make_tcm(
                tii,
                remove_cols=False,
                remove_rows=remove_rows,
                remove_sources=remove_sources,
                remove_ncs=remove_ncs,
            )
            tcmlist.append(tcm)
            latlist.append(np.array(model_lat))
            lonlist.append(np.array(model_lon))
        t3 = np.concatenate(tcmlist, axis=0)
        lat = np.concatenate(latlist, axis=0)
        lon = np.concatenate(lonlist, axis=0)
        self.latlist = np.array(latlist)
        self.lonlist = np.array(lonlist)
        if remove_cols:
            nmax = t3.shape[1]
            iremove = []
            # very last column is obs. So start with second to last column.
            for nnn in np.arange(nmax - 2, 0, -1):
                test = t3[:, nnn]
                if np.all(test == 0.0):
                    iremove.append(nnn)
                else:
                    break
            t3 = np.delete(t3, iremove, axis=1)
            self.tcm_columns = np.delete(self.tcm_columns, iremove, axis=0)
        self.tcm = t3
        self.tcm_lat = lat
        self.tcm_lon = lon
        self.latlist = np.array(latlist)
        self.lonlist = np.array(lonlist)
        return t3, lat, lon

    def remove_near_clear_sky(self, avg, window):
        """ """
        # InverseAsh class
        # this creates rolling average so nearby 0 pixels will have above zero values.
        avg2 = avg.rolling(x=window, center=True).max()
        avg2 = avg2.rolling(y=window, center=True).max()
        # areas above 0 in the smeared obs are set to True
        test1 = xr.where(avg2 > 0, True, False)
        # areas above 0 in the original obs are set to True
        test2 = xr.where(avg > 0, True, False)
        # areas from the original array set back to original value.
        # above 0 areas from smeared array are set to 0.
        # zero values are set to -1.
        test3 = xr.where(np.any([test1, test2], axis=0), avg, -1)
        # Returns xarray with
        # original above zero observations.
        # value of 0 in areas near to the above zero observations
        # value of -1 in areas far from the observations.
        return test3

    def make_tcm(
        self,
        tii,
        remove_cols=True,
        remove_rows=False,
        remove_sources=None,
        remove_ncs=0,
    ):
        # InverseAsh class
        """
        remove sources should be a list of  times / heights to remove
        along the ensemble dimension.
        Example.
        ensemble dimension names are generally "102119_2880" (MonthDayHour_BottomHeight)
        remove_sourecs = ['12880'] would remove all sources with _12880.
        """
        # remove rows means remove rows with observations of 0 or nan.

        # header line indicates release times for each column.
        # I think header can just be a dummy and it is added when writing to file.
        # one column for each release time/location.
        # one row for each measurement location.

        #         ens1  ens2 ens3 ens4 ens5 ens6 ens7 .... Obs
        # (x1,y1)
        # (x2,y1)
        # (x3,y1)
        # .
        # .
        # .
        # (x1,y2)

        # column will be from the ensemble dimension.
        # measurement is from the lat/lon dimension.

        # last column is the value of the observation.

        # number of rows corresponds to number of points where there is an observation
        # only include observations above 0.
        # or include where either observation OR model above 0.
        # or include only where model above 0.

        # right now only one time period.
        cdump = self.cdump_hash[tii]

        # remove some sources from consideration.
        if remove_sources:
            ekeep = cdump.ens.values
            for yyy in remove_sources:
                ekeep = [x for x in ekeep if yyy not in x]
            cdump = cdump.sel(ens=ekeep)

        avg = self.volcat_avg_hash[tii]

        s1 = avg.shape[0] * avg.shape[1]
        if remove_ncs > 0:
            avg = self.remove_near_clear_sky(avg, remove_ncs)

        cdump = cdump * self.concmult
        model = cdump.stack(pos=["y", "x"])
        model = model.transpose("pos", "ens")
        # some have nans? Find out why?
        model = model.fillna(0)
        # remove columns which have no contribution at all.
        if remove_cols:
            model = model.where(model > 0)
            model = model.dropna(dim="ens", how="all")

        model_lat = model.latitude.values.reshape(s1, 1)
        model_lon = model.longitude.values.reshape(s1, 1)
        columns = model.ens.values

        model = model.values

        volc = avg.values.reshape(s1, 1)
        volc_lat = avg.latitude.values.reshape(s1, 1)
        volc_lon = avg.longitude.values.reshape(s1, 1)

        volc = volc.flatten()

        if remove_ncs > 0:
            # remove only values that are 0.
            vpi = np.where(volc != 0)
            model = model[vpi]
            volc = volc[vpi]
            model_lat = model_lat[vpi]
            model_lon = model_lon[vpi]
            volc_lon = volc_lon[vpi]
            volc_lat = volc_lat[vpi]
            # set values that are less than 0 back to 0.
            volc = xr.where(volc < 0, 0, volc)

        if remove_rows:
            # only consider rows where observations are greater than 0.
            vpi = np.where(volc > 0)
            model = model[vpi]
            volc = volc[vpi]
            model_lat = model_lat[vpi]
            model_lon = model_lon[vpi]
            volc_lon = volc_lon[vpi]
            volc_lat = volc_lat[vpi]

        volc = volc.reshape(volc.shape[0], 1)

        tcm = np.concatenate([model, volc], axis=1)

        if not np.all(volc_lon == model_lon):
            print("WARNING, model and observed locations in tcm not matching")
        if not np.all(volc_lat == model_lat):
            print("WARNING, model and observed locations in tcm not matching")

        self.tcm = tcm
        self.tcm_lat = model_lat
        self.tcm_lon = model_lon
        # this contains the keys that can be matched in the sourcehash attribute.
        self.tcm_columns = columns

        return tcm, model_lat, model_lon, columns

    def make_tcm_columns(
        self,
        tii,
        remove_cols=True,
        remove_rows=False,
        remove_sources=None,
        remove_ncs=0,
    ):
        # just get the columns without creating the whole TCM.
        cdump = self.cdump_hash[tii]

        # remove some sources from consideration.
        if remove_sources:
            ekeep = cdump.ens.values
            for yyy in remove_sources:
                ekeep = [x for x in ekeep if yyy not in x]
            cdump = cdump.sel(ens=ekeep)

        # cdump = cdump * self.concmult
        model = cdump.stack(pos=["y", "x"])
        model = model.transpose("pos", "ens")
        # some have nans? Find out why?
        model = model.fillna(0)

        columns = model.ens.values
        # self.tcm_columns = columns
        return columns

    def plot_tcm(self):
        """ """
        # InverseAsh class
        cb = plt.pcolormesh(np.log10(self.tcm), cmap="tab20")
        plt.colorbar(cb)

    def make_fake(self):
        """
        not functional yet.
        """
        # InverseAsh class
        tcm_real = self.tcm.copy()
        nra = []
        nra.append([1, 0, 0, 0, 0, 0, 2])
        nra.append([0, 1, 0, 0, 0, 0, 3])
        nra.append([0, 0, 1, 0, 0, 0, 10])
        nra.append([0, 0, 0, 1, 0, 0, 6])
        nra.append([0, 0, 0, 0, 1, 0, 5])
        nra.append([0, 0, 0, 0, 0, 1, 4])
        self.tcm = np.array(nra)
        self.write_tcm("test_tcm.txt")
        self.tcm = tcm_real

    def write_tcm(self, name):
        """
        name : str
        """
        # InverseAsh class
        astr = ""
        sep = " "
        hstr = ""  # header string
        # print(self.tcm.shape)\
        # this is number of columns minus 1.
        # and needs to be input into Parameters.in.dat
        self.n_ctrl = self.tcm.shape[1] - 1
        # print('N_ctrl {}'.format(self.tcm.shape[1]-1))
        # print('output file {}'.format(name))
        for iii, line in enumerate(self.tcm):
            for jjj, val in enumerate(line):
                if iii == 0:
                    # write a dummy header line.
                    hstr += "43637.750" + sep
                if not np.isnan(val):
                    astr += "{:1.5e}".format(val)
                else:
                    astr += "{:1.4e}".format(0)
                astr += sep
            astr += "\n "
            if iii == 0:
                hstr += "\n"
        with open(name, "w") as fid:
            fid.write(hstr + astr)
        return hstr + astr

    def make_outdat(self, dfdat):
        # InverseAsh class
        """
        dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        Returns
        vals : tuple (date, height, emission mass)

        matches emissions from the out.dat file with
        the date and time of emission.
        uses the tcm_columns array which has the key
        and the sourehash dictionary which contains the information.
        """
        datelist = []
        htlist = []
        valra = []

        for val in zip(self.tcm_columns, dfdat[1]):
            shash = self.sourcehash[val[0]]
            datelist.append(shash["sdate"])
            htlist.append(shash["bottom"])
            valra.append(val[1])
        vals = list(zip(datelist, htlist, valra))
        return vals

    def make_outdat_df(self, dfdat, savename=None, part="basic"):
        # InverseAsh class
        """
        dfdat : pandas dataframe output by InverseOutDat class get_emis method.
                this is a list of tuples (source tag), value from emissions
        """
        vals = self.make_outdat(dfdat)
        vals = list(zip(*vals))
        ht = vals[1]
        time = vals[0]
        # emit = np.array(vals[2])/1.0e3/3600.0
        emit = np.array(vals[2])

        # this is for particle size. Used by the InverseAshPart class.
        if len(vals) == 4:
            psize = vals[3]
            data = zip(time, ht, psize, emit)
            if part == "index":
                iii = 0
                cols = [1, 2]
            elif part == "cols":
                iii = 1
                cols = [0, 2]
            colnames = ["date", "ht", "psize", "mass"]
        # this is for only height and time
        else:
            data = zip(time, ht, emit)
            iii = 1
            cols = 0
            colnames = ["date", "ht", "mass"]
        dfout = pd.DataFrame(data)
        if part == "condense" and len(vals) == 4:
            dfout.columns = colnames
            dfout = dfout.groupby(["date", "ht"]).sum()
            dfout = dfout.reset_index()
            dfout.columns = [0, 1, 2]
            iii = 1
            cols = 0

        if part == "basic":
            dfout.columns = colnames
            return dfout
        dfout = dfout.pivot(columns=cols, index=iii)
        if savename:
            print("saving  emissions to ", savename)
            dfout.to_csv(savename)
        return dfout

    def plot_outdat(
        self, vals, log=False, fignum=1, cmap="Blues", unit="kg/s", thresh=0
    ):
        # InverseAsh class
        """
        vals is output by make_outdat.
        """
        fig = plt.figure(fignum, figsize=(10, 5))
        vals = list(zip(*vals))
        sns.set()
        sns.set_style("whitegrid")
        # output in kg/s?/
        if unit == "kg/s":
            emit = np.array(vals[2]) / 1.0e3 / 3600.0
        elif unit == "kg/h":
            emit = np.array(vals[2]) / 1.0e3
        elif unit == "g/h":
            emit = np.array(vals[2]) / 1.0
        vpi = np.where(emit < thresh)
        emit[vpi] = 0
        ht = np.array(vals[1]) / 1e3
        if not log:
            cb = plt.scatter(vals[0], ht, c=emit, s=100, cmap=cmap, marker="s")
        else:

            cb = plt.scatter(
                vals[0], ht, c=np.log10(emit), s=100, cmap=cmap, marker="s"
            )
            # cb = plt.pcolormesh(vals[0],ht,emit,cmap=cmap)
        cbar = plt.colorbar(cb)
        cbar.ax.set_ylabel(unit)
        fig.autofmt_xdate()

    def plot_outdat_profile(
        self, dfdat, fignum=1, unit="kg", ax=None, clr="--ko", alpha=1, lw=1
    ):
        """ """
        # InverseAsh class
        if not ax:
            sns.set()
            sns.set_style("whitegrid")
            fig = plt.figure(fignum, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
        df = self.make_outdat_df(dfdat, part="cols")
        ts = df.sum(axis=1)

        sns.set()
        # the one represents one hour time period.
        yval = ts.values * 1 / 1e12
        ax.plot(
            yval, ts.index.values / 1000.0, clr, alpha=alpha, linewidth=lw, label="one"
        )
        ax.set_xlabel("Tg of mass emitted", fontsize=15)
        ax.set_ylabel("Height (km)", fontsize=15)
        totalmass = yval.sum()
        # print('total {} Tg'.format(totalmass))
        return ax, df

    def plot_outdat_ts(
        self, dfdat, log=False, fignum=1, unit="kg/s", ax=None, clr="--ko"
    ):
        # InverseAsh class
        """
        plots time series of MER. summed along column.
        dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        """
        if not ax:
            sns.set()
            sns.set_style("whitegrid")
            fig = plt.figure(fignum, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
        df = self.make_outdat_df(dfdat, part="index")
        sns.set()
        ts = df.sum()
        if unit == "kg/s":
            yval = ts.values / 3.6e6
        elif unit == "g/h":
            yval = ts.values
        # print('plotting outdat ts', yval, type(yval))
        # print([x[1] for x in ts.index.values])
        ax.plot([x[1] for x in ts.index.values], yval, clr, marker="o", label="one")
        # fig.autofmt_xdate()
        ax.set_ylabel("MER {}".format(unit), fontsize=15)
        handles, labels = ax.get_legend_handles_labels()
        return ax, df

    def plot_out2dat_times(self, df2, cmap="viridis"):
        """
        not working
        """
        # InverseAsh class
        return -1

    def plot_out2dat_scatter(self, tiilist, df2, vloc, cmap="Blues"):
        """ """
        # InverseAsh class
        if isinstance(tiilist, int):
            tiilist = [tiilist]
        sns.set()
        ppp = 0
        # tii = self.time_index(daterange[0])
        modelall = df2["model"].values
        nnn = 0
        for tii in tiilist:
            fig = plt.figure(figsize=[10, 5])
            ax1 = fig.add_subplot(1, 2, 1)
            ax2 = fig.add_subplot(1, 2, 2)
            volcat = self.volcat_avg_hash[tii]
            shape = volcat.shape
            lon = self.lonlist[tii - 1]
            lat = self.latlist[tii - 1]
            model = modelall[nnn : nnn + len(lon)]
            volcat = self.volcat_avg_hash[tii]
            r2 = volcat.where(volcat > 0)
            cb = ax1.scatter(lon, lat, c=model, cmap=cmap, s=10, marker="o")
            cb2 = ax2.scatter(
                volcat.longitude,
                volcat.latitude,
                c=r2.values,
                s=10,
                cmap=cmap,
                marker="o",
            )
            nnn = len(lon)
            if vloc:
                ax1.plot(vloc[0], vloc[1], "y^")
                ax2.plot(vloc[0], vloc[1], "y^")
            plt.colorbar(cb, ax=ax1)
            plt.colorbar(cb2, ax=ax2)

    def plot_out2dat(self, tiilist, df2, cmap="viridis", vloc=None, ptype="pcolormesh"):
        # InverseAsh class
        """
        tiilist needs to be same order as for tcm.
        df2 is from InverseOutDat get_conc method.
        doesn't work when only later time periods are used!
        """
        if isinstance(tiilist, int):
            tiilist = [tiilist]
        sns.set()
        ppp = 0
        # tii = self.time_index(daterange[0])
        for tii in tiilist:
            fig = plt.figure(figsize=[10, 5])
            ax1 = fig.add_subplot(1, 2, 1)
            ax2 = fig.add_subplot(1, 2, 2)
            volcat = self.volcat_avg_hash[tii]
            shape = volcat.shape
            model = df2["model"].values
            temp = model[ppp : ppp + shape[0] * shape[1]]
            model = temp.reshape(shape[0], shape[1])
            vpi = np.where(model < 0.01)
            model[vpi] = np.nan
            ppp = ppp + shape[0] * shape[1]
            r2 = volcat.where(volcat > 0)
            m_max = np.nanmax(model)
            v_max = np.nanmax(r2.values)
            m_min = np.nanmin(model)
            v_min = np.nanmin(r2.values)
            p_min = np.nanmin([m_min, v_min])
            p_max = np.nanmax([m_max, v_max])
            norm = mpl.colors.Normalize(vmin=p_min, vmax=p_max)
            print(np.nanmax(model), np.nanmax(r2.values))
            if ptype == "pcolormesh":
                cb = ax1.pcolormesh(
                    volcat.longitude,
                    volcat.latitude,
                    model,
                    norm=norm,
                    cmap=cmap,
                    shading="nearest",
                )
                cb2 = ax2.pcolormesh(
                    volcat.longitude,
                    volcat.latitude,
                    r2.values,
                    norm=norm,
                    cmap=cmap,
                    shading="nearest",
                )
            # cb = ax1.scatter(volcat.longitude, volcat.latitude,c=np.log10(model),cmap=cmap,s=50,marker='s')
            else:
                # lon = self.lonlist[tii-1][0:14]
                # lat = self.latlist[tii-1][0:21]
                # xv, yv = np.meshgrid(lon,lat)
                lon = self.lonlist[tii - 1][0:14]
                lat = self.latlist[tii - 1][0:21]
                xv, yv = np.meshgrid(lon, lat)
                print("vshape", volcat.longitude.shape)
                print(xv.shape, yv.shape)
                print("model shape", model.shape)
                cb = ax1.scatter(xv, yv, c=model, cmap=cmap, s=10, marker="o")
                cb2 = ax2.scatter(
                    volcat.longitude,
                    volcat.latitude,
                    c=r2.values,
                    s=10,
                    cmap=cmap,
                    marker="o",
                )
            plt.colorbar(cb, ax=ax1)
            # cb2 = ax2.scatter(volcat.longitude, volcat.latitude, c=r2.values,s=10,cmap=cmap,marker='o')
            plt.colorbar(cb2, ax=ax2)
            if vloc:
                ax1.plot(vloc[0], vloc[1], "y^")
                ax2.plot(vloc[0], vloc[1], "y^")
            plt.show()
        return volcat

    def compare_plotsA(self, daterange, zii=None, tii=None, levels=None):
        # InverseAsh class
        """
        must input either daterange or tii.

        if zii is None then sum along ensemble dimension showing coverage of all HYSPLIT runs.
        For the inversion runs, the release from different heights is shown by the ens dimension.

        """

        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(1, 1, 1)
        if tii is None:
            tii = self.time_index(daterange[0])
        print("tii", tii)
        cdump = self.concmult * self.cdump_hash[tii]
        volcat = self.volcat_avg_hash[tii]
        if not zii:
            csum = cdump.sum(dim="ens")
        else:
            csum = cdump.isel(ens=zii)
            print(cdump.ens.values[zii])
            # print(csum.ens.values)
            # print(self.sourcehash[str(csum.ens.values)])
        # print(cdump.time)
        # print(csum.coords)
        # volcat.plot.pcolormesh(x='longitude',y='latitude',levels=levels,ax=ax1)
        # cdump.sum(dim='ens').plot.contour(x='longitude',y='latitude',ax=ax2)
        try:
            # plt.pcolormesh(csum.longitude, csum.latitude, np.log10(csum),cmap='Reds',shading='nearest')
            cbm = ax1.pcolormesh(
                csum.x, csum.y, np.log10(csum), cmap="Reds", shading="nearest"
            )
        except:
            print("FAILED max value", np.max(csum))
            print("------------------")
        # plt.pcolormesh(csum.longitude, csum.latitude, csum,cmap='Reds',shading='nearest')
        # cb= csum.plot.pcolormesh(x='longitude',y='latitude',cmap='viridis',ax=ax2)
        cb = ax1.pcolormesh(
            volcat.x,
            volcat.y,
            np.log10(volcat),
            cmap="Blues",
            shading="nearest",
            alpha=0.5,
        )
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=np.log10(volcat),s=2,cmap='Blues')
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=volcat.values,s=2,cmap='viridis',levels=levels)
        # vals = np.log10(volcat)
        # cb = plt.contour(volcat.x, volcat.y, vals.fillna(0),cmap='viridis',levels=[0,1,10,100])
        plt.colorbar(cbm)
        plt.tight_layout()
        # return csum.copy()
        return csum.copy()

    def get_norm(self, model, r2, logscale=False):
        # InverseAsh class
        """ """
        if isinstance(r2, np.ndarray):
            rval = r2
        else:
            rval = r2.values

        m_max = np.nanmax(model)
        v_max = np.nanmax(rval)
        m_min = np.nanmin(model)
        v_min = np.nanmin(rval)
        p_min = np.nanmin([m_min, v_min])
        p_max = np.nanmax([m_max, v_max])
        lower_color = 0.8
        if logscale:
            norm = mpl.colors.LogNorm(vmin=lower_color * p_min, vmax=p_max)
        else:
            norm = mpl.colors.Normalize(vmin=p_min, vmax=p_max)
        return norm

    # def generate_pairs(self):
    #    for tii in self.volcat_avg_hash.keys():
    #        volcat = self.volcat_avg_hash[iii]
    #        cdump = self.cdump_hash[iii]*self.concmult
    #        yield volcat, cdump

    def set_bias_correction(self, slope, intercept, dfcdf=pd.DataFrame()):
        # InverseAsh class
        """ """
        self.slope = slope
        self.intercept = intercept
        self.dfcdf = dfcdf

    def manage_cdump(
        self,
        tii,
        coarsen=None,
        slope=None,
        intercept=None,
        dfcdf=pd.DataFrame(),
        cii=None,
        coarsen_max=None,
    ):
        # InverseAsh class
        """
        tii : datetime.datetime OR integer
        coarsen :  int : if > 0 then use the xarray coarsen method to coarsen the data using mean.
        coarsen_max :  int : if > 0 then use the xarray coarsen method to coarsen the data using max.

        bias correction.
        For bias correction, first will try to use the dfcdf DataFrame if input.
        If it is empty will check to see if the dfcdf attribute is not empty.
        If the dfcdf attribute is not empty then it will use that.
        If the dfcdf attribute and dfcdf input are both empty then will
        use the slope and intercept inputs. If those are both None then no bias correction is applied.

        cii   : datetime.datetime. used with dfcdf
        dfcdf : pandas DataFrame with information for bias correction.
        slope : float : if not None used for bias correction if  dfcdf is empty
        intecept : float : if not None used for bias correction if  dfcdf is empty
        """

        # make sure that tii is datime and iii is integer - index for tii.
        if isinstance(tii, int):
            iii = tii
            tii = self.get_time(iii)
        elif isinstance(tii, datetime.datetime):
            iii = self.time_index(tii)

        cdump = self.cdump_hash[iii] * self.concmult

        if not cii:
            cii = tii

        # cii = tii-datetime.timedelta(hours=2)
        # if pandas dataframe available use that.
        if dfcdf.empty and not self.dfcdf.empty:
            dfcdf = self.dfcdf
        if not dfcdf.empty:
            # create xarray with slope and intercept with coordinate of ens.
            temp = dfcdf[dfcdf["time"] == cii]
            temp = temp[["slope", "intercept", "ens"]]
            temp = temp.set_index("ens")
            xrt = temp.to_xarray()
            # apply slope and intercept correction
            cdump_prev = cdump.copy()
            # print('cdump A', cdump.ens.values)
            cdump = cdump * (1 - xrt.slope)
            cdump = cdump - xrt.intercept
            f2 = xr.where(cdump > 0, cdump, 0)
            # print('xrt', xrt.ens.values)
            # print('f2', f2.ens.values)
            # print('cdump', cdump.ens.values)
            # print('prev', cdump_prev.ens.values)
            f2 = xr.where(cdump_prev < 0.0001, 0, f2)
            cdump = f2
        # otherwise use slope and intercept if those are specified.
        else:
            if not slope and self.slope:
                slope = self.slope
            cdump_prev = cdump
            if slope:
                # print('multiplying cdump by {}'.format(1-slope))
                cdump = cdump * (1 - slope)
            if not intercept and self.intercept:
                intercept = self.intercept
            if intercept:
                # print('shifting cdump by {}'.format(-1*intercept))
                cdump = cdump - intercept
                # if intercept is positive then remove negative values.
                f2 = xr.where(cdump > 0, cdump, 0)
                # if intercept is negative then do not add to values that were previously 0.
                f2 = xr.where(cdump_prev < 0.0001, 0, f2)
                cdump = f2

        # coarsen comes last.
        # find the mean value
        if coarsen:
            cdump = cdump.coarsen(x=coarsen, boundary="trim").mean()
            cdump = cdump.coarsen(y=coarsen, boundary="trim").mean()
        # use maximum value in scene
        elif coarsen_max:
            cdump = cdump.fillna(0)
            cdump = cdump.coarsen(x=coarsen_max, boundary="trim").max()
            cdump = cdump.coarsen(y=coarsen_max, boundary="trim").max()

        return cdump

    def get_pair(
        self,
        tii,
        coarsen=None,
        slope=None,
        intercept=None,
        dfcdf=pd.DataFrame(),
        cii=None,
        coarsen_max=None,
    ):
        # InverseAsh class
        """ """
        if isinstance(tii, int):
            iii = tii
        elif isinstance(tii, datetime.datetime):
            iii = self.time_index(tii)
        volcat = self.volcat_avg_hash[iii]
        cdump = self.manage_cdump(
            tii, coarsen, slope, intercept, dfcdf, cii, coarsen_max
        )
        if coarsen:
            volcat = volcat.coarsen(x=coarsen, boundary="trim").mean()
            volcat = volcat.coarsen(y=coarsen, boundary="trim").mean()
        elif coarsen_max:
            volcat = volcat.coarsen(x=coarsen_max, boundary="trim").max()
            volcat = volcat.coarsen(y=coarsen_max, boundary="trim").max()
        return volcat, cdump

    def match_volcat(self, forecast):
        # InverseAsh class
        """ """
        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        volcat = self.volcat_avg_hash[tii]
        return volcat

    def compare_forecast_dist(self, forecast, thresh=None):
        # InverseAsh class
        """ """
        from utilhysplit.evaluation import statmain

        volcat = self.match_volcat(forecast)
        evals = forecast.values
        vvals = volcat.values
        exval, eyval = statmain.nancdf(evals.flatten(), thresh)
        vxval, vyval = statmain.nancdf(vvals.flatten(), thresh)
        fig = plt.figure(1)
        ax = fig.add_subplot(1, 1, 1)
        ax.step(exval, eyval, "-r", label="HYSPLIT")
        ax.step(vxval, vyval, "-b", label="Volcat")
        ks1, ks2 = statmain.kstestnan(evals.flatten(), vvals.flatten(), thresh)
        try:
            print(
                "Kolmogorov-Smirnov Parameter {} {}".format(
                    np.max(np.abs(ks1)), np.max(np.abs(ks2))
                )
            )
        except:
            return 1
        return np.max(np.abs(ks1))

    def remove_nans(self, data, thresh=None):
        # InverseAsh class
        """ """
        # remove nans
        data2 = data[~np.isnan(data)]
        if thresh:
            data2 = data2[data2 > thresh]
            # apply threshold.
            # vpi = data2 < thresh
            # data2[vpi] = np.nan
            # data2 = data2[~np.isnan(data2)]
        return data2

    def compare(self, forecast, thresh=None):
        # InverseAsh class
        """ """
        volcat = self.match_volcat(forecast)
        evals = self.remove_nans(forecast.values.flatten(), thresh)
        vvals = self.remove_nans(volcat.values.flatten(), thresh)
        # number of pixles above threshold in each.
        forecast_parea = len(evals)
        volcat_parea = len(vvals)
        return forecast_parea, volcat_parea

    def compare_forecast_alt(
        self, forecast, cmap="Blues", ptype="pcolormesh", vloc=None
    ):
        # InverseAsh class
        """ """
        # forecast should be an xarray in mass loading format with no time dimension.
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(figsize=[15, 5])
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3)

        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        print("tii", tii)
        volcat = self.volcat_avg_hash[tii]
        evals = forecast.values
        vpi = evals < 0.001
        evals[vpi] = np.nan
        norm = self.get_norm(volcat, forecast)

        # evals = np.log10(evals)
        # vvals = np.log10(volcat.values)
        vvals = volcat.values
        vpi = vvals < 0.001
        vvals[vpi] = np.nan
        # clevels = [0.2,2,5,10]
        clevels = [0.2, 10]
        if ptype == "pcolormesh":
            cb = ax1.pcolormesh(
                volcat.longitude,
                volcat.latitude,
                vvals,
                norm=norm,
                cmap=cmap,
                shading="nearest",
            )
            cb2 = ax2.pcolormesh(
                forecast.longitude,
                forecast.latitude,
                evals,
                norm=norm,
                cmap=cmap,
                shading="nearest",
            )
            # cb2 = ax3.pcolormesh(forecast.longitude, forecast.latitude,evals,norm=norm, cmap=cmap,shading='nearest')
            ax3.contourf(
                volcat.longitude,
                volcat.latitude,
                volcat.values,
                levels=clevels,
                cmap="Reds",
            )
            ax3.contour(
                forecast.longitude,
                forecast.latitude,
                evals,
                levels=clevels,
                cmap="viridis",
            )
        plt.title(time.strftime("%Y %m/%d %H:%M UTC"))
        plt.colorbar(cb, ax=ax1)
        plt.colorbar(cb2, ax=ax2)
        ylim = ax1.get_ylim()
        ax2.set_ylim(ylim)
        ax3.set_ylim(ylim)
        xlim = ax1.get_xlim()
        ax2.set_xlim(xlim)
        ax3.set_xlim(xlim)
        if vloc:
            ax1.plot(vloc[0], vloc[1], "y^")
            ax2.plot(vloc[0], vloc[1], "y^")
            ax3.plot(vloc[0], vloc[1], "y^")
        return fig

    def compare_forecast(self, forecast, cmap="Blues", ptype="pcolormesh", vloc=None):
        # InverseAsh class
        """
        forecast should be an xarray in mass loading format with no time dimension.
        """
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(figsize=[15, 5])
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3)

        time = pd.to_datetime(forecast.time.values)
        tii = self.time_index(time)
        print("tii", tii)
        volcat = self.volcat_avg_hash[tii]
        evals = forecast.values
        vpi = evals < 0.001
        evals[vpi] = np.nan
        norm = self.get_norm(volcat, forecast)

        # evals = np.log10(evals)
        # vvals = np.log10(volcat.values)
        vvals = volcat.values
        vpi = vvals < 0.001
        vvals[vpi] = np.nan
        # clevels = [0.2,2,5,10]
        clevels = [0.2, 10]
        if ptype == "pcolormesh":
            cb = ax1.pcolormesh(
                volcat.longitude,
                volcat.latitude,
                vvals,
                norm=norm,
                cmap=cmap,
                shading="nearest",
            )
            cb2 = ax2.pcolormesh(
                forecast.longitude,
                forecast.latitude,
                evals,
                norm=norm,
                cmap=cmap,
                shading="nearest",
            )
            # cb2 = ax3.pcolormesh(forecast.longitude, forecast.latitude,evals,norm=norm, cmap=cmap,shading='nearest')
            ax3.contourf(
                volcat.longitude,
                volcat.latitude,
                volcat.values,
                levels=clevels,
                cmap="Reds",
            )
            ax3.contour(
                forecast.longitude,
                forecast.latitude,
                evals,
                levels=clevels,
                cmap="viridis",
            )
        plt.title(time.strftime("%Y %m/%d %H:%M UTC"))
        plt.colorbar(cb, ax=ax1)
        plt.colorbar(cb2, ax=ax2)
        ylim = ax1.get_ylim()
        ax2.set_ylim(ylim)
        ax3.set_ylim(ylim)
        xlim = ax1.get_xlim()
        ax2.set_xlim(xlim)
        ax3.set_xlim(xlim)
        if vloc:
            ax1.plot(vloc[0], vloc[1], "y^")
            ax2.plot(vloc[0], vloc[1], "y^")
            ax3.plot(vloc[0], vloc[1], "y^")
        return fig

    def compare_plots(self, daterange, levels=None):
        # InverseAsh class
        """ """
        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        tii = self.time_index(daterange[0])
        cdump = self.concmult * self.cdump_hash[tii]
        volcat = self.volcat_avg_hash[tii]
        csum = cdump.sum(dim="ens")
        volcat.plot.pcolormesh(x="longitude", y="latitude", levels=levels, ax=ax1)
        # cdump.sum(dim='ens').plot.contour(x='longitude',y='latitude',ax=ax2)
        # plt.pcolormesh(csum.longitude, csum.latitude, np.log10(csum),cmap='Reds')
        # plt.pcolormesh(csum.longitude, csum.latitude, csum,cmap='Reds',shading='nearest')
        cb = csum.plot.pcolormesh(x="longitude", y="latitude", cmap="viridis", ax=ax2)
        # cb = plt.pcolormesh(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues',levels=levels)
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=np.log10(volcat),s=2,cmap='Blues')
        # cb = plt.scatter(volcat.longitude, volcat.latitude, c=volcat.values,s=2,cmap='viridis',levels=levels)
        # cb = plt.contour(volcat.longitude, volcat.latitude, np.log10(volcat),cmap='Blues')
        # plt.colorbar(cb)
        plt.tight_layout()
        return ax1, ax2

    def prepare_times(self, daterangelist, das=None):
        # InverseAsh class
        """ """
        for dates in daterangelist:
            subdas = []
            for sd in das:
                if (
                    pd.to_datetime(sd.time.values) < dates[1]
                    and pd.to_datetime(sd.time.values) >= dates[0]
                ):
                    subdas.append(sd)
            print("preparing", dates)
            self.prepare_one_time(dates, das=subdas)

    def add_volcat_hash(self, volcat_hash):
        # InverseAsh class
        """ """
        # allow to update from another instance of the class.
        self.volcat_avg_hash.update(volcat_hash)

    def add_cdump_hash(self, cdump_hash):
        # InverseAsh class
        """ """
        self.cdump_hash.update(cdump_hash)

    def change_sampling_time(self, value):
        if value not in ["start", "end"]:
            print("time description not recognized {}".format(value))
            print("value should be {} or {}".format("start", "end"))
        else:
            self.cdump.attrs.update({"time description": value})

    def get_sampling_time(self, default="start"):
        if "time description" in self.cdump.attrs.keys():
            if "start" in str.lower(self.cdump.attrs["time description"]):
                return "start"
            if "end" in str.lower(self.cdump.attrs["time description"]):
                return "end"
        else:
            return default

    def mass_load_modified(self):
        # logic for only getting mass loading in certain levels.

        # have an array of top heights.
        #
        return -1

    def coarsen(self, num=3):
        if not self.coarsen:
            self.original_cdump_hash = self.cdump_hash.copy()
            self.original_volcat_avg_hash = self.volcat_avg_hash.copy()
        self.coarsen = num
        for key in self.cdump_hash.keys():
            c2 = self.original_cdump_hash[key].coarsen(x=num).mean()
            c2 = c2.coarsen(y=num).mean()
            v2 = self.original_volcat_avg_hash[key].coarsen(x=num).mean()
            v2 = v2.coarsen(y=num).mean()
            self.cdump_hash[key] = c2
            self.volcat_hash[key] = v2

    def prepare_one_time(self, daterange, htoptions=0, zvals=None, verbose=False):
        """
        daterange : list of datetime objects
        htoptions : boolean : capability not fully implemented.
        zvals : list of vertical levels (in meters) to use to calculate mass loading.


        Populates the folling dictionaries:
        self.cdump_hash[tii] = cdump_a
        self.volcat_avg_hash[tii] = avra
        self.volcat_ht_hash[tii] = maxvhra
        """

        # currently must coarsen all data at the same time.
        if self.coarsen:
            print(
                "Warning: Adding new data after some data has already been coarsened."
            )
        vdir = self.vdir
        vid = self.vid
        # key for hashes is determined from times in cdump file.
        # check whether cdump time is start or end of sampling time.
        if self.get_sampling_time() == "start":
            # use first date in daterange list
            model_tii = 0
        elif self.get_sampling_time() == "end":
            # use second date in daterange list
            model_tii = 1
        # if st == "end":
        #    model_tii = 1
        # else:
        #    model_tii = 0
        tii = self.time_index(daterange[model_tii])
        if not isinstance(tii, int):
            print("No time found for {}".format(daterange[0]))
            return None, None
        if isinstance(zvals, np.ndarray):
            zvals = list(zvals)
        cdump_a = hysplit.hysp_massload(
            self.cdump.sel(time=daterange[model_tii]), zvals=zvals
        )
        if htoptions == 1:
            cdump_b = self.cdump.sel(time=daterange[model_tii])

        vset = self.get_volcat(daterange, verbose)
        buf = 5
        # clip the volcat array before aligning.
        try:
            vra = vset.ash_mass_loading
        except:
            return None, None

        # also get the volcat observed height.
        try:
            vhra = vset.ash_cloud_height
        except:
            return None, None
        vset.close()
        a1, a2, b1, b2 = self.clip(vra.sum(dim="time"), buf=buf)
        vra = vra.transpose("time", "y", "x")
        # vra has dimensions of time, y, x
        vra = vra[:, a1:a2, b1:b2]

        vhra = vhra.transpose("time", "y", "x")
        vhra = vhra[:, a1:a2, b1:b2]

        # clip the cdump array before aligning.
        if "ens" in cdump_a.coords:
            cdump_a = cdump_a.transpose("ens", "y", "x")
            if htoptions == 1:
                cdump_b = cdump_b.transpose("ens", "z,", "y", "x")
            dummy = cdump_a.sum(dim="ens")
        else:
            cdump_a = cdump_a.transpose("y", "x")
            if htoptions == 1:
                cdump_b = cdump_b.transpose("z", "y", "x")
            dummy = cdump_a

        try:
            a1, a2, b1, b2 = self.clip(dummy, buf=5)
        except:
            print("dummy cannot be clipped", dummy)
        if "ens" in cdump_a.coords:
            cdump_a = cdump_a[:, a1:a2, b1:b2]
            if htoptions == 1:
                cdump_b = cdump_b[:, :, a1:a2, b1:b2]
        else:
            cdump_a = cdump_a[a1:a2, b1:b2]
            if htoptions == 1:
                cdump_b = cdump_b[:, a1:a2, b1:b2]

        # align the grids.
        tol = 0.01
        if not hgu.compare_grids(cdump_a, vra, verbose=False, tolerance=tol):
            if not hgu.compare_grids_compat(cdump_a, vra):
                avra = vra.fillna(0).mean(dim="time")
                self.volcat_avg_hash[tii] = avra
                print("prepare_one_time: grids are not compatible")
                return False
            if verbose:
                print("grids are compatible. changing grids")
            vatt = hgu.find_grid_specs(vra)
            dlat = vatt["Latitude Spacing"]
            dlon = vatt["Longitude Spacing"]
            crnr = vatt["llcrnr longitude"]
            newatt = hgu.create_global_grid(dlat, dlon, crnr)
            cdump_a = hgu.change_grid(cdump_a, newatt)
            vhra = hgu.change_grid(vhra, newatt)
            # print(hgu.find_grid_specs(cdump_a))
            # print(vatt)
        if hgu.compare_grids(cdump_a, vra, verbose=True, tolerance=tol):
            dummy, vhra = hgu.align_grids(cdump_a, vhra, tolerance=tol)
            cdump_a, vra = hgu.align_grids(cdump_a, vra, tolerance=tol)
        else:
            print("----------")
            hgu.compare_grids(cdump_a, vra, verbose=True)
            print("prepare_one_time: grids cannot be aligned")
            print(daterange)
            # print(vset)
            print("----------")
            return cdump_a, vra

        # take the average of the volcat data.
        avra = vra.fillna(0).mean(dim="time")
        maxvhra = vhra.fillna(0).max(dim="time")

        self.cdump_hash[tii] = cdump_a
        self.volcat_avg_hash[tii] = avra
        self.volcat_ht_hash[tii] = maxvhra
        return cdump_a, avra

    def compare_time_ave(self, daterange):
        """
        creates plots illustrating how time averaging affects
        volcat data.
        """
        sns.set()
        vset = self.get_volcat(daterange)
        vset = vset.ash_mass_loading
        for time in vset.time.values:
            temp = vset.sel(time=time)
            a1, a2, b1, b2 = self.clip(temp.fillna(0))
            temp = temp[a1:a2, b1:b2]
            temp.plot.pcolormesh(x="longitude", y="latitude")
            plt.show()
        avra = vset.fillna(0).mean(dim="time")
        avra2 = vset.mean(dim="time")
        a1, a2, b1, b2 = self.clip(avra)
        avra = avra[a1:a2, b1:b2]
        avra2 = avra2[a1:a2, b1:b2]
        print("Average with nans set to 0")
        avra.plot.pcolormesh(x="longitude", y="latitude")
        plt.show()
        print("Average of values above 0")
        avra2.plot.pcolormesh(x="longitude", y="latitude")
        plt.show()
        diff = avra2.fillna(0) - avra
        diff.plot.pcolormesh(x="longitude", y="latitude")


def make_control(
    efile,
    tdir="./",
    sname="CONTROL.default",
    wdir="./",
    suffix="emit",
    metstr="gfs0p25",
    lat=0,
    lon=0,
    forecast_directory=None,
    archive_directory=None,
    duration=None,
):

    """
    duration : if set will over-ride run duration in CONTROL.default.
    """
    # set up met data finder.
    metfilefinder = MetFileFinder(metstr.lower())
    metfilefinder.set_forecast_directory(forecast_directory)
    metfilefinder.set_archives_directory(archive_directory)
    if "gefs" in metstr.lower():
        metfilefinder.set_ens_member("." + suffix)

    control = hcontrol.HycsControl(fname=sname, working_directory=tdir)
    control.read()
    # number of locations need to be written.
    if duration:
        control.run_duration = duration
    nlocs = efile.cycle_list[0].nrecs / len(efile.splist)
    nspecies = len(efile.sphash.keys())
    stime = efile.sdatelist[0]
    done = False
    while not done:
        if nlocs == control.nlocs:
            done = True
        else:
            control.add_dummy_location()
    # set to 0 to use emit times file.
    for species in control.species:
        species.rate = 0
        species.duration = 0
    for grid in control.concgrids:
        grid.outfile = "cdump.{}".format(suffix)
        grid.outdir = wdir + "/"
        grid.centerlat = float(np.floor(lat))
        grid.centerlon = float(np.floor(lon))
    # update metfiles.
    control.remove_metfile(rall=True)
    metfiles = metfilefinder.find(stime, int(control.run_duration), hours=24)
    for mfile in metfiles:
        control.add_metfile(mfile[0], mfile[1])

    # simulation start date same as emit-times
    control.rename("CONTROL.{}".format(suffix), wdir)
    control.add_sdate(stime)
    control.write()
    return control


def make_setup(
    tdir="./", sname="SETUP.default", wdir="./", efilename="emit.txt", suffix="emit"
):
    print("making setup")
    setup = hcontrol.NameList(sname, tdir)
    setup.read()
    if suffix == "emit":
        setup.add("efile", efilename)
    else:
        setup.add("efile", "{}_{}".format(suffix, efilename))
    setup.add("poutf", "PARDUMP.{}".format(suffix))
    setup.rename("SETUP.{}".format(suffix), wdir)
    setup.write()
    return setup


def construct_efile(
    vals,
    vlat,
    vlon,
    area=1,
    emis_threshold=50,
    vres=1000,
    name="emit.txt",
    date_cutoff=None,
    phash=None,
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
    from utilhysplit import emitimes

    # print("construct_efile function")
    efile = emitimes.EmiTimes(species=phash)
    # print("splist", efile.splist)
    # efile.set_species(phash)
    duration = "0100"
    cycle_list = []

    addso4 = False
    if "pSO4" in phash.keys():
        addso4 = True

    for iii, value in enumerate(vals):
        time = value[0]
        if date_cutoff:
            if time >= date_cutoff:
                break
        height = value[1]
        emis = value[2]
        if len(value) > 3:
            part = phash[value[3]]
        else:
            part = 1
        # print("construct_efile function:", value, part)
        if emis < emis_threshold:
            emis = 0
        newcycle = False
        if time not in cycle_list:
            newcycle = True
        if newcycle:
            efile.add_cycle(time, duration=1)
            cycle_list.append(time)
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
    phash=None,
):
    efile = construct_efile(
        vals, vlat, vlon, area, emis_threshold, vres, name, date_cutoff, phash
    )
    efile.write_new(name)
    return efile


def plot_outdat_profile_psize_function(
    dfdat, log=False, fignum=1, unit="kg/s", ax=None, cmap="viridis"
):
    # plots time series of MER. summed along column.
    # dfdat : pandas dataframe output by make_outdat_df function with
    #        part='basic'
    if not ax:
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(fignum, figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
    if isinstance(dfdat, str):
        df = pd.read_csv(dfdat, index_col=0, header=1)
        df = df.dropna()
    else:
        df = dfdat
    sns.set()
    plist = df.psize.unique()
    cm = colormaker.ColorMaker(cmap, len(plist), ctype="hex", transparency=None)
    cmlist = cm()
    for iii, psize in enumerate(plist):
        dfp = df[df.psize == psize]
        dfp = dfp.drop("psize", axis=1)
        # dfp = dfp.pivot(index='date',columns='ht')
        dfp = dfp.pivot(index="date", columns="ht")
        try:
            dfp = dfp.mass
        except:
            pass
        ts = dfp.sum()
        xval = ts.index.values * 1 / 1.0e3
        yvals = ts.values / 1.0e12
        ax.plot(yvals, xval, "#" + cmlist[iii], label=psize)
    ax.set_xlabel("Tg of mass emitted", fontsize=15)
    ax.set_ylabel("Height (km)", fontsize=15)
    # ax.set_ylabel('Mass {}'.format(unit),fontsize=15)
    return ax, ts


def plot_outdat_ts_psize_function(
    dfdat, log=False, fignum=1, unit="kg/s", ax=None, cmap="viridis"
):
    # plots time series of MER. summed along column.
    # dfdat : pandas dataframe output by make_outdat_df function with
    #        part='basic'

    if not ax:
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(fignum, figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
    if isinstance(dfdat, str):
        df = pd.read_csv(dfdat, index_col=0, header=1)
        df = df.dropna()
    else:
        df = dfdat
    sns.set()
    plist = df.psize.unique()
    cm = colormaker.ColorMaker(cmap, len(plist), ctype="hex", transparency=None)
    cmlist = cm()
    for iii, psize in enumerate(plist):
        dfp = df[df.psize == psize]
        dfp = dfp.drop("psize", axis=1)
        # dfp = dfp.pivot(index='date',columns='ht')
        dfp = dfp.pivot(index="ht", columns="date")
        try:
            dfp = dfp.mass
        except:
            pass
        ts = dfp.sum()
        if unit == "kg/s":
            yval = ts.values / 3.6e6
        elif unit == "g/h":
            yval = ts.values
        print(type(ts.index.values[0]), ts.index.values[0])
        xval = [pd.to_datetime(x) for x in ts.index.values]
        # ax.plot([x[0] for x in ts.index.values], yval, clr)
        print(cmlist[iii])
        ax.plot(xval, yval, "#" + cmlist[iii], label=psize)
        # fig.autofmt_xdate()
        ax.set_ylabel("MER {}".format(unit), fontsize=15)
    return ax, ts


def plot_outdat_ts_function(
    dfdat,
    log=False,
    fignum=1,
    unit="kg/s",
    ax=None,
    clr="--ko",
    label=None,
    alpha=1,
    lw=1,
    marker=None,
):
    # plots time series of MER. summed along column.
    # dfdat : pandas dataframe output by InverseOutDat class get_emis method.
    if not ax:
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(fignum, figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
    if isinstance(dfdat, str):
        df = pd.read_csv(dfdat, index_col=0, header=1)
        df = df.dropna()
    else:
        df = dfdat
    sns.set()
    ts = df.sum()
    if unit == "kg/s":
        yval = ts.values / 3.6e6
    elif unit == "g/h":
        yval = ts.values
    xval = [pd.to_datetime(x) for x in ts.index.values]
    # ax.plot([x[0] for x in ts.index.values], yval, clr)
    ax.plot(
        xval, yval, label=label, color=clr, linestyle="-", marker=marker, alpha=alpha, linewidth=lw
    )
    # fig.autofmt_xdate()
    ax.set_ylabel("MER {}".format(unit), fontsize=15)
    return ax, df


def plot_outdat_profile_function(
    dfdat, fignum=1, unit="km", ax=None, clr="k", label=None, alpha=1, lw=1, marker="o"
):
    if not ax:
        sns.set()
        sns.set_style("whitegrid")
        fig = plt.figure(fignum, figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)

    if isinstance(dfdat, str):
        df = pd.read_csv(dfdat, index_col=0, header=1)
        df = df.dropna()
    else:
        df = dfdat

    ts = df.sum(axis=1)
    sns.set()
    xval = ts.values * 1 / 1e12
    try:
        yval = ts.index.values / 1000.0
    except:
        yval = list(map(float, list(ts.index.values)))
        yval = np.array(yval) / 1000.0
    if unit == "FL":
        yval = yval * 3280.84 / 100.0
    ax.plot(
        xval,
        yval,
        color=clr,
        marker=marker,
        label=label,
        alpha=alpha,
        linewidth=lw,
        linestyle="-",
    )
    ax.set_xlabel("Tg of mass emitted", fontsize=15)
    if unit == "FL":
        ax.set_ylabel("Height (FL)", fontsize=15)
    else:
        ax.set_ylabel("Height (km)", fontsize=15)
    totalmass = xval.sum()
    # print('total {} Tg'.format(totalmass))
    return ax, totalmass


class InverseSO2(InverseAsh):
    def prepare_one_time(
        self, daterange, das=None, htoptions=0, st="start", zvals=None
    ):
        return -1


class InverseAshPart(InverseAsh):
    def get_cdump(self, tdirlist, fnamelist, verbose, dummy):
        # vdir,
        # vid,
        # configdir='./',configfile=None,
        # ensdim='ens',
        # verbose=False):
        psizelist = ["p060", "p200", "p500"]
        spnum = np.arange(1, len(psizelist) + 1)
        xrlist = []

        for fname in zip(tdirlist, fnamelist):
            print("appending ", fname)
            temp = xr.open_dataset(os.path.join(fname[0], fname[1]))
            temp = temp[list(temp.data_vars.keys())[0]]
            temp = temp.rename("pall")
            xrlist.append(temp)
        xrlist = xrlist
        self.cdump = stack_psizes(xrlist, psizelist)
        self.phash = self.add_phash()
        for xrl in xrlist:
            xrl.close()

    def add_phash(self, phash=None):
        """
        phash : dictionary

        """
        if not phash:
            ptemp = list(set(self.cdump.ens.psize.values))
            ptemp.sort()
            iii = np.arange(1, len(ptemp) + 1)
            self.phash = dict(zip(ptemp, iii))
        return self.phash

    def make_outdat(self, dfdat):
        """
        makeoutdat for InverseAshPart class.
        dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        Returns
        vals : tuple (date, height, emission mass, particle size)
        """
        # matches emissions from the out.dat file with
        # the date and time of emission.
        # uses the tcm_columns array which has the key
        # and the sourehash dictionary which contains the information.
        datelist = []
        htlist = []
        valra = []
        psizera = []
        # if not isinstance(self.tcm_columns, (list, np.ndarray)):
        #   self.tcm_columns = self.make_tcm_columns()
        for val in zip(self.tcm_columns, dfdat[1]):
            shash = self.sourcehash[val[0][0]]
            datelist.append(shash["sdate"])
            htlist.append(shash["bottom"])
            valra.append(val[1])
            psizera.append(val[0][1])
        vals = list(zip(datelist, htlist, valra, psizera))
        return vals

    def plot_outdat_ts(self, dfdat, log=False, fignum=1, unit="kg/s", ax=None):
        # InverseAshPart class

        # plots time series of MER. summed along column.
        # dfdat : pandas dataframe output by InverseOutDat class get_emis method.
        if not ax:
            sns.set()
            sns.set_style("whitegrid")
            fig = plt.figure(fignum, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
        df = self.make_outdat_df(dfdat, part="basic")
        ax, ts = plot_outdat_ts_psize_function(df)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        return ax, ts
        # sns.set()
        # ts = df.sum()
        # if unit == 'kg/s':
        #   yval = ts.values/3.6e6
        ##jjelif unit == 'g/h':
        #   yval = ts.values
        # ax.plot([x[1] for x in ts.index.values], yval, clr)
        # fig.autofmt_xdate()
        # ax.set_ylabel('MER {}'.format(unit),fontsize=15)
        # return ax, df

    def plot_outdat_profile(self, dfdat, fignum=1, unit="kg", ax=None):
        # InverseAshPart class
        if not ax:
            sns.set()
            sns.set_style("whitegrid")
            fig = plt.figure(fignum, figsize=(10, 5))
            ax = fig.add_subplot(1, 1, 1)
        df = self.make_outdat_df(dfdat, part="basic")
        ax, ts = plot_outdat_profile_psize_function(df)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        return ax, ts
