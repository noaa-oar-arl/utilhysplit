#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# OutputDisersion
#
# -----------------------------------------------------------------------------
# Creates netcdf file xrfilejobid.nc from cdump files.
# TODO add capabability to create the AWIPS2 formatted netcdf files.
# this is currently in ashbase.py
# -----------------------------------------------------------------------------


# from abc mport ABC, abstractmethod
import logging
import os

from ashapp.ashruninterface import ModelOutputInterface
from ashapp import utils
from ashapp.ashnetcdf import HYSPLITAshNetcdf


# from runhandler import ProcessList
from utilvolc.volcMER import HT2unit

logger = logging.getLogger(__name__)
utils.setup_logger()




class OutputDispersion(ModelOutputInterface):

    ilist = [('fraction_of_fine_ash','req'),
             ('eflag','req'),
             ('Use_Mastin_eq','req'),  # whether to use mastin eq. to calculate mer.
             ('start_date','req') 
            ]


    def __init__(self, inp, filelist):
        """
        filelist: list of tuples (cdumpname, metfile, sourceid)
        """


        self.JOBID = inp["jobid"]
        self.inp = inp
        self._inputlist = filelist
        self._outputlist = self.set_outputlist()

        self._ncfile = HYSPLITAshNetcdf(self.cdumpxraname())


    @property
    def ncfile(self):
        return self._ncfile

    @property
    def inputlist(self):
        return self._inputlist

    @inputlist.setter
    def inputlist(self, inlist):
        self._inputlist = inlist

    # no setter. read only
    @property
    def outputlist(self):
        return self._outputlist

    def cdumpxraname(self):
        return os.path.join(self.inp['WORK_DIR'], "xrfile.{}.nc".format(self.JOBID)) 
    def set_outputlist(self):
        outputlist = []
        outputlist.append(self.cdumpxraname())
        return outputlist

    def postprocess(self):
        self.get_cdump_xra()
        self._ncfile.write_with_compression(overwrite=False)
        # combine all cdump files into an xarray and write a netcdf file.
        return -1

    def check(self):
        for filename in self._outputlist:
            if not os.path.isfile(filename):
                return False
        return True

    ## additional methods ----------------------------------------------------
    def get_cdump_xra(self):
        """
        reads from file if it exists.
        Otherwise create from list of cdump names.
        updates the multiplication factor if necessary.
        """
        if self._ncfile.empty():
            blist = []
            blist = [x for x in self.inputlist if "cdump" in x[0]]
            century = 100 * (int(self.inp["start_date"].year / 100))
            species = None
            ainp = self.inp.copy()
            ainp["mult"] = 1
            if 'polygon' in ainp.keys():
               ainp['polygon'] = str(ainp['polygon'])
            self._ncfile.make_cdump_xra(blist, century, species=species, inp=ainp)
        mult = self.get_conc_multiplier()
        change = self._ncfile.changemult(mult)
        return self._ncfile.cxra

    def get_conc_multiplier(self):
        """
        factor to convert concentration in unit mass /m3 to mg/m3.
        Model must have been run with 1 unit mass/h released.
        """
        if self.inp['Use_Mastin_eq']:
            height = (self.inp['top'] - self.inp['bottom']) / 1000.0
            # HT2unit outputs in grams. Then multiply by 1e3 to get mg
            conc_multiplier = 1e3 * HT2unit(height) * self.get_ash_reduction()
        else:
            conc_multiplier = self.get_ash_reduction()
        return conc_multiplier

    def get_ash_reduction(self):
        """
        """
        # 09/07/2023 (amc) don't allow multiplication by 0.
        eflag = float(self.inp["eflag"])
        if eflag==0: eflag = 1
        M63 = self.inp['fraction_of_fine_ash']  # fraction of  fine ash
        # 07/01/2023 decided to chang eflag use.
        # it was done this way to mimic the ash reduction VAAC currently use.
        #conc_multiplier = 10 ** (-1 * eflag)
        conc_multiplier = eflag * M63
        return M63 * conc_multiplier


