# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# ashruninterface.py - run HYSPLIT model on web and create plots
#
#
# -----------------------------------------------------------------------------
# Provides class interfaces
#
# -----------------------------------------------------------------------------

from abc import ABC, abstractmethod



class ModelRunCollection(ABC):
         
    @property
    @abstractmethod
    def inp(self):
      """
      dictionary with user inputs
      """
      pass

    @property
    @abstractmethod
    def processhandler(self):
      """
      ProcessList class
      """
      pass

    #@property
    #@abstractmethod
    #def model_list(self):
    #  """
    #  list of  ModelRunInterface objects
    #  """
    #  pass

    @abstractmethod
    def run_model(self):
      """
      Run the model
      """
      pass


class ModelRunInterface(ABC):
    # ModelRun consists of the following objects
    # CONTROL file
    # SETUP file
    # metfilefinder
    # filelocator/filename composer
    # dictionary of inputs
    # 

    @abstractmethod
    def compose_control(stage,rtype):
      """
      Create the control file
      """
      pass

    @abstractmethod
    def compose_setup(stage):
      """
      Create the SETUP.CFG file
      """
      pass

    @abstractmethod
    def create_run_command(self):
      """
      creates command str
      """
      pass

    @abstractmethod
    def run_model(self):
      """
      complete setup for model to run.
      """
      pass

    # Properties.

    @property
    @abstractmethod
    def control(self):
      """
      control file hcontrol.HycsControl object
      """
      pass

    @property
    @abstractmethod
    def filelist(self):
      """
      list of names of  files created
      """
      pass

    @property
    @abstractmethod
    def filelocator(self):
      """
      FileLocator object which creates filenames
      """
      pass

    @property
    @abstractmethod
    def inp(self):
      """
      dictionary with user inputs
      """
      pass

    @property
    @abstractmethod
    def metfilefinder(self):
      """
      MetFileFinder object
      """
      pass

    @property
    @abstractmethod
    def setup(self):
      """
      setup file hcontrol.Namelist object
      """
      pass

    @property
    @abstractmethod
    def status(self):
      """
      returns tuple of (str, list)
      str indicates current status
          values can be INITIALIZED, FAILED, or COMPLETE
      list is a list of strings with history of
          status and other actions.

      """
      pass

class ModelOutputInterface(ABC):

    @abstractmethod
    def postprocess(self):
        """
        post-processing on model output
        """
        pass

    @abstractmethod
    def check(self):
        """
        check to see all model output exists
        """
        pass

    @property
    @abstractmethod
    def inputlist(self):
      """
      list input files
      """
      pass

    @property
    @abstractmethod
    def outputlist(self):
      """
      list of output files
      """
      pass

class MainRunInterface(ABC):
    #def __init__(sele, JOBID):
    # Ash Run consists of the following components
    # 0. inputs
    # 1. model run(s)
    # 2. model output
    # 3. graphics
    # Methods are
    # add_api_info
    # update_run_status
    # handle_crash

    @abstractmethod
    def add_api_info(self, apistr, urlstr, headerstr):
        pass

    @abstractmethod
    def update_run_status(self, jobId, status):
        pass

    @abstractmethod
    def handle_crash(self, stage=0):
        pass

    @abstractmethod
    def debug_message(self):
        pass

    @abstractmethod
    def doit(self):
        """ 
        main method
        """
        pass

    @property
    @abstractmethod
    def inp(self):
        """
        Dictionar with inputs
        """
        pass

    @property
    @abstractmethod
    def JOBID(self):
        """
        string with job identification 
        """
        pass

    @property
    @abstractmethod
    def modelgraphics(self):
        """
        ModelOutputInterface class
        """
        pass

    @property
    @abstractmethod
    def modeloutput(self):
        """
        ModelOutputInterface class
        """
        pass

    @property
    @abstractmethod
    def modelrun(self):
        """
        ModelRunInterface class
        """
        pass


