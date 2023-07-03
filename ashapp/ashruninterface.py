# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# ashruninterface.py - run HYSPLIT model  and postprocess output.
#
# -----------------------------------------------------------------------------
# Provides class interfaces for running model and postprocessing model output.

# First the model runs must be setup and run.

# ModelCollectionInterface
#    Interface for more than one run
#    EnsembleDispersion - uses the RunDispersion class.
#    CollectInverse - uses RunDispersion class
#    CollectEmitTimes - uses RunEmitTimes

# ModelRunInterface
#    Interface for creating the model runs
#    RunDispersion - used for regular dispersion runs
#                  - used for unit source runs for inversion
#    RunEmitTimes - used for data insertion 
#                 - TODO also use for running results of inversion.
#    TODO - RunTrajectory,
#           

# Then the output must be processed.

# ModelOutputInterface
#     Interface for processing model output.
#     OutputDispersion

# -----------------------------------------------------------------------------
# The OutputDispersion class is utilized by all in order to
# put information from all model runs into a netcdf file using hysplit.combine_dataset.

# Simple dispersion runs
#   MainDispersion(MainRunInterface)
#     RunDispersion(ModelRunInterface)
#     OutputDispersion(ModelOutputInterface)
#     TODO - graphics

# Ensemble dispersion runs
#   MainEnsemble(MainDispersion)
#     EnsembleDispersion(ModelCollectionInterface)
#     OutputDispersion(ModelOutputInterface)
#     TODO - grapics

# Data insertion runs
#   MainEmitTimes(MainDispersion)
#     CollectEmitTimes(ModelCollectionInterface)
#     OutputDispersion(ModelOutputInterface)
#     TODO -graphics  

# Inversion runs using source term
#     MainEmitTimes(MainDispersion)

# Inversion unit source runs
#   MainInverse(MainDispersion)
#     CollectInverse(ModelCollectionInterface)
#     OutputDispersion(ModelOutputInterface)
#     TODO - graphics

# Inversion unit source runs with GEFS
#     MainGEFSInverse(MainInverse)
#     CollectInverse(ModelCollectionInterface)
#     OutputDispersion(ModelOutputInterface)
#     TODO - graphics

# TODO
#     forward trajectory runs
#     ensemble trajectory runs
#     backward trajectory runs


#-------------
# ash_main.py
#
# ensembledispersion.py
# collectemittimes.py
# collectinverse.py

# graphicsdispersion.py (TODO- needs work)

# rundispersion.py
# runemittimes.py

# maindispersion.py
# mainensemble.py

# outputdispersion.py
# 
# HELPER files
# ashnetcdf.py - helper for the OutputDispersionClass




from abc import ABC, abstractmethod



class ModelRunInterface(ABC):
    """ 
    #  ModelRun consists of the following objects
    # CONTROL file
    # SETUP file
    # metfilefinder
    # filelocator/filename composer
    # dictionary of inputs
    """
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
    def run(self,overwrite):
      """
      call the run command
      """
      pass


    @abstractmethod
    def run_model(self,overwrite):
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
    """
     The Main Run Interface 
      consists of the following properties.
     0. inputs
     1. model run(s) 
     2. model output
        a. list of files to be available for download.
     3. model output graphics
        a. list of files to be available for download. 

     Methods are
        add_api_info
        update_run_status
        handle_crash
    """
    # where does the after run check go?
 
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




class ModelCollectionInterface(ABC):


    @property
    @abstractmethod
    def filehash(self):
      """
      More comprehensive list of files
      """
      pass

    @property
    @abstractmethod
    def filelist(self):
      """
      List of main outputs (cdump files)
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
    def status(self):
      """
      status
      """
      pass

    #@property
    #@abstractmethod
    #def processhandler(self):
    #  """
    #  ProcessList class
    #  """
    #  pass

    #@property
    #@abstractmethod
    #def model_list(self):
    #  """
    #  list of  ModelRunInterface objects
    #  """
    #  pass

    @abstractmethod
    def setup(self):
      """
      Run the model
      """
      pass

    @abstractmethod
    def run(self,overwrite):
      """
      Run the model
      """
      pass


