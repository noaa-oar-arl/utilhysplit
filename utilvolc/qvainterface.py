# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# qvainterface.py - run HYSPLIT model on web and create plots
#
#
# -----------------------------------------------------------------------------
# Provides class interfaces
# DFInterface - interface for managing a dataframe associated with a class
# StatusInterface - 
# -----------------------------------------------------------------------------
from abc import ABC, abstractmethod


class AshParamsInterface(ABC):

    @property
    @abstractmethod
    def gridspace(self):
        """
        grid spacing. used for HYSPLIT runs concentration grid
        as well as write_parallax_corrected for regridding volcat data.
        """
        pass




class DFInterface(ABC):
    #TODO -  may want to add data types for the columns as well.
 
    @property
    @abstractmethod
    def edf(self):
        """
        pandas dataframe
        """
        pass

    @property
    @abstractmethod
    def required_columns(self):
        """
        list of required column names
        """
        pass

    @abstractmethod
    def are_required_columns_present(self):
        pass

    @abstractmethod
    def add_df(self):
        """
        method to add another dataframe.
        """
        pass

    @abstractmethod
    def read(self):
        """
        reads a previously saved csv file and adds it to dataframe
        """
        pass

    @abstractmethod
    def save(self):
        """
        saves the dataframe as a csv file
        """
        pass


class EventInterface(ABC):

    @property
    @abstractmethod
    def edf(self):
        pass

    @abstractmethod
    def set_dir(self):
        """
        set directories for 
        downloading volcat netcdf files.
        writing emittimes files
        writing parallax corrected files
        doing inversion runs                
        """
        pass

    @abstractmethod
    def get_dir(self):
        pass

    @abstractmethod
    def set_volcano_name(self):
        pass
 
    @abstractmethod
    def download(self):
        """
        download volcat netcdf files
        """
        pass
 
    @abstractmethod
    def write_emit(self):
        """
        write emit-times for data insertion
        """
        pass

    @abstractmethod
    def write_parallax_corrected(self):
        """
        write parallax corrected files
        """
        pass

    @abstractmethod
    def InvInfo(self):
        """
        Basic information such as estimated start time, plume height etc.
        """
        pass

    


    #------------PLOTS---------------------------------------------
    # plots(pste,pc=False,levels,vlist,central_longitude)
    # plots_with_vaas(self,vaas,pstep,pc=True)
    # compare_pc(pstep,daterange,fid,central_longitude,vlist)
    # vplot
    # boxplot

    @abstractmethod
    def compare_pc(self):
        """
        compare parallax corrected to non-parallax corrected
        """
        pass

    @property
    @abstractmethod
    def df(self):
        """
        Eventdf class
        """
        pass

    @property
    @abstractmethod
    def vloc(self):
        """
        list [latitude, longitude]
        """
        pass

    @property
    @abstractmethod
    def ndir(self):
        """
        str
        """
        pass

    @property
    @abstractmethod
    def events(self):
        """
        list of xarray DataSets with volcat data
        """
        pass
 
    @property
    @abstractmethod
    def pcevents(self):
        """
        list of xarray DataSets with volcat data
        """
        pass
