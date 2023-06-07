import logging
from ashapp.ashruninterface import ModelRunCollection
from ashapp.rundispersion import RunDispersion
from utilhysplit.metfiles import gefs_suffix_list


logger = logging.getLogger(__name__)

class EnsembleDispersion(ModelRunCollection):

    def __init__(self,inp,jobid):
        self.JOBID=jobid

        self._ilist = ['meteorologicalData','forecastDirectory','archivesDirectory',
                 'WORK_DIR','HYSPLIT_DIR','jobname','durationOfSimulation','latitude',
                 'longitude','bottom','top','emissionHours','rate','area','start_date',
                 'samplingIntervalHours','jobid']

        self._inp = {}
        self.inp = inp
        
        self._filehash = {}
        self._filelist = []

        self._status = {'MAIN':'INITIALIZED'}

    @property
    def filehash(self):
        return self._filehash 

    @property
    def filelist(self):
        return self._filelist

    @property
    def inp(self):
        return self._inp

    @inp.setter
    def inp(self,inp):
        self._inp.update(inp)
        complete = True
        for iii in self._ilist:
            if iii not in self._inp.keys(): 
               logger.warning('Input does not contain {}'.format(iii))
               complete=False
        if 'jobid' in self._inp.keys():
            self.JOBID = self._inp['jobid'] 
        if complete: logger.info('Input contains all fields')

    @property
    def processhandler(self):
        return self._processhandler

    @property
    def status(self):
        return self._status

    def setup_runs(self):
        inp = self.inp.copy()
        command_list = []
        for suffix in gefs_suffix_list():
            inp['jobid'] = '{}_{}'.format(self.JOBID,suffix)
            run = RunDispersion(inp)
            run.metfilefinder.set_ens_member("." + suffix)
            command = run.run_model(overwrite=False)
            self._status[suffix] = run.status
            if 'FAILED' in run.status[0] or 'COMPLETE' in run.status[0]:
               logger.warning(run.status)
               continue
            if command: 
               command_list.append(command)
            self._filehash.update(run.filehash)
            self._filelist.extend(run.filelist)
        return command_list

    def run_model(self):
        command_list = self.setup_runs()



 
