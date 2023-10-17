# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# inversioninterface.py - inverse modeling
#
#
# -----------------------------------------------------------------------------
# Provides class interfaces

# InversionInterface
# TCMInterface
# PairedDataInterface
# RunInversionInterface


from abc import ABC, abstractmethod


class InversionInterface(ABC):

   @property
   @abstractmethod
   def taglist(self):
       """
       all the unit source runs from the model
       """  
       pass

   @abstractmethod
   def close_arrays(self):
        pass    


class TCMInterface(ABC):

   @property
   @abstractmethod
   def columns(self):
       """
       The column names
       """  
       pass


   @abstractmethod
   def write(self,name):
       pass 

   @abstractmethod
   def make_tcm_mult(self,tiilist,
                     remove_cols,
                     remove_rows,
                     remove_sources,
                     remove_ncs):
       pass 

   @abstractmethod
   def make_tcm(self,tii,
                 remove_cols,
                 remove_rows,
                 remove_sources,
                 remove_ncs):
       pass 

class PairedDataInterface(ABC):
   """
   Paired data has modeloutput and observations.
   """
 
   @property
   @abstractmethod
   def inp(self):
       """
       dictionary with configuration information
       """  
       pass

   @property
   @abstractmethod
   def modeloutput(self):
       """
       all the unit source runs from the model
       """  
       pass

   @property
   @abstractmethod
   def obs(self):
       """
       Observations
       """  
       pass

   @abstractmethod
   def close_arrays(self):
       pass

   @abstractmethod
   def add_model(self):
       """
       add netcdf representation of cdump file.
       """
       pass

   @abstractmethod
   def prepare_one_time(self):
       """
       matches the observations and model data for one time period.
       """
       pass

   @abstractmethod
   def coarsen(self):
       """
       takes the otuput from prepare_one_time and creates
       a spatial average
       """
       pass


class RunInversionInterface(ABC):
   """
   RunInversion has  paired data, a tcm,
   """
 
   @property
   @abstractmethod
   def directories(self):
      pass
 
   @property
   @abstractmethod
   def inp(self):
       """
       dictionary with configuration information
       """  
       pass
 
   @property
   @abstractmethod
   def paired_data(self):
       """
       all the unit source runs from the model
       """  
       pass

   @property
   @abstractmethod
   def tcm(self):
       """
       Transfer coefficient matrix
       """  
       pass

   @abstractmethod
   def copy(self):
       pass

   @abstractmethod
   def print_summary(self):
       pass


