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

from monetio.models import hytraj


logger = logging.getLogger(__name__)
utils.setup_logger()


class OutputTrajectory(ModelOutputInterface):

    #ilist = [('starte_date','req')]
    ilist = ['jobid']

    def __init__(self, inp, filelist):
        """
        filelist: list of tuples (cdumpname, metfile, sourceid)
        """
        self.JOBID = inp["jobid"]
        self.inp = inp
        self._inputlist = filelist
        self._outputlist = self.set_outputlist()

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

    def set_outputlist(self):
        outputlist = [self.what_is_outputfilename()]
        return outputlist

    def postprocess(self):
        flist = []
        taglist = []
        for tdump in self._inputlist:
            flist.append(tdump[0])
            taglist.append(tdump[1])
            pass 
        # trajectory  dataframe
        tdf =  hytraj.combine_dataset(flist,taglist,renumber=True)
        tdf.to_csv(self.what_is_outputfilename())
        return tdf

    def what_is_outputfilename(self):
        return "tdump.{}.csv".format(self.JOBID)

    def check(self):
        for filename in self._outputlist:
            if not os.path.isfile(filename):
                return False
        return True

    ## additional methods ----------------------------------------------------

