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
        return "xrfile.{}.nc".format(self.JOBID)

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
            # cdumpname = cdumpname[0]
            # logger.info('Creating xra for cdump files {}'.format(cdumpname))
            # source_tag = "Line to {:1.0f} km".format(self.inp["top"] / 1000.0)
            # met_tag = self.inp["meteorologicalData"]
            # blist = [(cdumpname, source_tag, met_tag)]
            century = 100 * (int(self.inp["start_date"].year / 100))
            species = None
            ainp = self.inp.copy()
            ainp["mult"] = 1
            # print(type(ainp["generatingPostscript"]))
            # ainp goes into the attributes for the netcdf file.
            self._ncfile.make_cdump_xra(blist, century, species=species, inp=ainp)
        mult = self.get_conc_multiplier()
        change = self._ncfile.changemult(mult)
        return self._ncfile.cxra

    # no longer needed?
    #def write_cdump_xra(self):
        # convert to mg /m3
   #     mult = self.get_conc_multiplier()
   #     change = self._ncfile.changemult(mult)
   #     if change:
   #         ainp = self.inp
   #         ainp["mult"] = mult
   #         Helper.remove(os.path.join(self.inp["WORK_DIR"], fname))
   #         self._ncfile.assign_attrs(ainp)
   #         self._ncfile.write_with_compression()
#
#        if not self._ncfile.empty():
#            logger.info(
#                "make_awips_netcdf: cxra empty. cannot create awips\
#                         files"
#            )
#            return []

    # TODO modify this.
    def make_awips_netcdf(self):
        #import cdump2netcdf
        # TODO update this!
        #ghash = {}
        #ghash["source_latitude"] =  self.inp["latitude"]
        #ghash["source_longitude"] = self.inp["longitude"]
        #ghash["source_name"] =      self.inp["VolcanoName"]
        #ghash["emission_start"] =   self.inp["start_date"]
        #ghash["emission_duration_hours"] = self.inp["emissionHours"]
        # in mg
        #mult = self.get_conc_multiplier(self.inp['top'],self.inp['bottom'])
        #mer = mult / 1e6 / 3600  # kg released in a second.
        #ghash["MER"] = mer
        #ghash["MER_unit"] = "kg/s"
        #logger.debug("MULT value for awips {:3e}".format(mult))
        #awipsname = self.filelocator.get_awips_filename(stage=0)
        #c2n = cdump2netcdf.Cdump2Awips(
        #    self._cxra, awipsname, munit="mg", jobid=self.JOBID, globalhash=ghash
        #)
        #awips_files = c2n.create_all_files()
        awips_files = []
        # returns list of awips files that were created.
        return awips_files

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
