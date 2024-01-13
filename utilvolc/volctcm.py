#import datetime
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from utilvolc.runhelper import Helper

from utilvolc.inversioninterface import TCMInterface
from utilvolc.utiltcm import ParametersIn, InverseDat
from utilvolc import utiltcm

logger = logging.getLogger(__name__)

# classes
#   TCM

# Functions
#   make_outdat
#   make_outdat_df
#   remove_near_clear_sky


# 2023 Dec 04 (amc) added emtpy method to TCM.
# 2023 Dec 04 (amc) made n_ctrl a property
# 2023 Dec 04 (amc) added get_output method

class TCM(TCMInterface):
    """
    This class only handles the tcm and running the inversion algorithm.
    """

    def __init__(self, tag='TCM'):
        self._columns = None
        self._tcm = np.array([[]])

        self.tcm_name = 'TCM_sum.csv'

        self.emissions =  None #InvEstimatedEmissions class
        self.eval= None        #InverseOut2Dat class

        self.tag = tag

        self.latlist = None
        self.lonlist = None
        self.tcm_lat = None
        self.tcm_lon = None
        self._history = ["initialized"]

    @property
    def history(self):
        return self._history

    def empty(self):
        if self.n_ctrl==-1:
           return True
        else:
           return False

    @property
    def tcm(self):
        """
        """
        return self._tcm

    @property
    def columns(self):
        """
        These are the column headers
        """
        return self._columns
 
    @property
    def n_ctrl(self):
        return(self.tcm.shape[1] - 1)

    @columns.setter
    def columns(self, columns):
        self._columns = columns

    def make_tcm_mult(
        self,
        paired_data_list,
        concmult=1,
        remove_cols=True,
        remove_rows=True,
        remove_sources=None,
        remove_ncs=0,
    ):
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
        for pdata in paired_data_list:
            model = pdata[0]
            obs = pdata[1]
            # print(self.cdump.time.values[tii])
            tcm, model_lat, model_lon, columns = self.make_tcm(
                model, obs,
                concmult = concmult,
                remove_cols=False,
                remove_rows=remove_rows,
                remove_sources=remove_sources,
                remove_ncs=remove_ncs,
            )
            tcmlist.append(tcm)
            #print(type(model_lon),model_lon.shape, model_lat.shape)
            latlist.append(np.array(model_lat))
            lonlist.append(np.array(model_lon))
        t3 = np.concatenate(tcmlist, axis=0)
        lat = np.concatenate(latlist, axis=0)
        lon = np.concatenate(lonlist, axis=0)
        #print(latlist, type(latlist))
        #self.latlist = np.array(latlist)
        self.latlist = lat
        self.lonlist = lon  
        #self.lonlist = np.array(lonlist)
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
            if len(self.columns) > 0:
                self.columns = np.delete(self.columns, iremove, axis=0)
        self._tcm = t3
        self.tcm_lat = lat
        self.tcm_lon = lon
        #self.latlist = np.array(latlist)
        #self.lonlist = np.array(lonlist)
        print('HERE C')
        return t3, lat, lon

    #def make_outdat(self, sourcehash):
    #    dfdat = self.get_emis()
    #    return make_outdat(sourcehash, self.columns, dfdat)

    def plot(self):
        """ """
        cb = plt.pcolormesh(np.log10(self.tcm), cmap="tab20")
        plt.colorbar(cb)


    def write(self, name='TCM_sum.csv'):
        """
        name : str
        """
        astr = ""
        sep = " "
        hstr = ""  # header string
        # print(self.tcm.shape)\
        # this is number of columns minus 1.
        # and needs to be input into Parameters.in.dat
        columns = list(self.columns.copy())
        columns.append('obs')
        for iii, line in enumerate(self.tcm):
            for jjj, val in enumerate(line):
                if iii == 0:
                    # write the column header
                    if jjj <= len(columns)-1:
                        # Must be floats or integers in header row
                        #hstr += columns[jjj] + sep
                        hstr += str(jjj) + sep
                    else:
                        hstr += "-999" + sep
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
        self.tcm_name = name
        print('writing tcm to {}'.format(self.tcm_name))
        self._history.append('written')
        return hstr + astr

    def make_tcm(
        self,
        cdump,
        obsavg,
        concmult,
        remove_cols=True,
        remove_rows=False,
        remove_sources=None,
        remove_ncs=0,
    ):
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
        # cdump = self.cdump_hash[tii]

        # remove some sources from consideration.
        if remove_sources:
            ekeep = cdump.ens.values
            for yyy in remove_sources:
                ekeep = [x for x in ekeep if yyy not in x]
            cdump = cdump.sel(ens=ekeep)

        avg = obsavg

        s1 = avg.shape[0] * avg.shape[1]
        if remove_ncs > 0:
            avg = remove_near_clear_sky(avg, remove_ncs)

        cdump = cdump * concmult
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

        volc = volc.reshape(volc.shape[0], 1)

        tcm = np.concatenate([model, volc], axis=1)

        if not np.all(volc_lon == model_lon):
            print("WARNING, model and observed locations in tcm not matching")
        if not np.all(volc_lat == model_lat):
            print("WARNING, model and observed locations in tcm not matching")

        self._tcm = tcm
        self.tcm_lat = model_lat
        self.tcm_lon = model_lon
        # this contains the keys that can be matched in the sourcehash attribute.
        self.tcm_columns = columns
        return tcm, model_lat, model_lon, columns

    def make_tcm_names(self):
        #InverseAshEns 
        out_name1 = "out.dat"
        out_name2 = "out2.dat"
        name1 = "{}_{}".format(self.tag, out_name1)
        name2 = "{}_{}".format(self.tag, out_name2)
        return name1, name2

    #def get_output(self,subdir):
    #    new_name1, new_name2 = self.make_tcm_names()
    #    return InverseDat(wdir=subdir, fname = new_name1, fname2 = new_name2)

    def run(self,execdir,subdir):
        """
        """
        out_name1 = "out.dat"
        out_name2 = "out2.dat"
        inp_name = "TCM_sum.csv"
        cmd = os.path.join(execdir, "new_lbfgsb.x")
        new_name1, new_name2 = self.make_tcm_names()
        #for iii, tcm in enumerate(self.tcm_names):

        print("run_tcm tag", self.tag)
        os.chdir(subdir)
        print("working in ", os.getcwd())
        params = ParametersIn(os.path.join(execdir, "Parameters_in.dat.original"))
        params.change_and_write(
            os.path.join(subdir, "Parameters_in.dat"), self.n_ctrl
        )


        Helper.remove(inp_name)

        # remove any existing output files.
        Helper.remove(out_name1)
        Helper.remove(out_name2)

        # run
        Helper.copy(self.tcm_name, inp_name)
        # Helper.execute_with_shell(cmd)
        Helper.execute(cmd)
            # move output files to new names

        Helper.move(out_name1, new_name1)
        Helper.move(out_name2, new_name2)

        with open('fort.188','r') as fid:
             temp = fid.readlines()
        self.cost = [x for x in temp if 'Cost' in x]
        self.gradient = [x for x in temp if 'Gradient' in x] 
        self.result = temp[-10:]
        

        Helper.move("fort.188", "fort.188.{}".format(self.tag))
        try:
            self.create_emissions(subdir)
        except Exception as eee:
            logger.warning('could not get output {}'.format(eee))

    def create_emissions(self,subdir):
        new_name1, new_name2 = self.make_tcm_names()
        fname = os.path.join(subdir,new_name1)
        self.emissions= utiltcm.InvEstimatedEmissions(fname,columns=self.columns)
        self.emissions.read(subdir)


    def create_outdat(self,subdir):
        new_name1, new_name2 = self.make_tcm_names()
        fname = os.path.join(subdir,new_name2)
        self.out2dat = utiltcm.InverseOut2Dat(fname)
        self.out2dat.read(subdir)


    #def get_emissions_df(self,sourcehash):
    #    self.emissions.sourcehash = sourcehash
    #    return self.emissions.make_emissions()
        
 



# TODO also in volcinverse as a method.
def make_outdat_df(vals, savename=None, part="basic"):
    """
    dfdat : pandas dataframe output by InvEstimatedEmissions class get_emis method.
            this is a list of tuples (source tag), value from emissions
    """
    #vals = make_outdat(tcm_columns, dfdat)
    #vals = list(zip(*vals))
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




def remove_near_clear_sky(avg, window):
    """
    this creates rolling average so nearby 0 pixels will have above zero values.
    """
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
