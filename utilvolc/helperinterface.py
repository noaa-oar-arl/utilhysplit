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


class FileNameInterface(ABC):


    def __str__(self):
        """
        returns filename string
        """

    def parse(self,fname):
        """
        given a string representing filename, parse it for information
        which is stored in attributes.
        """
        pass

    def make_filename(self):
        """
        given some inputs create filename
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


