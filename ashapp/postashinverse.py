import datetime
import numpy as np
import os
import utilvolc.ash_eval as ae
import utilvolc.ash_inverse as ai
import glob
from monetio.models import hysplit
import matplotlib.pyplot as plt
from utilvolc.runhelper import Helper


def get_gefs_flist(fname):
    from utilhysplit.metfiles import gefs_suffix_list
    fname.replace('.nc','')
    ens = gefs_suffix_list()
    return ['{}_{}.nc'.format(fname,x) for x in ens]


class RunInversion:

    def __init__(self, tdirlist, fnamelist, vdir, vid, configdir, configfile):
        """
        tdirlist : list of str
        fnamelist : list of str
        vdir : str
        vid  : str
        configdir : str
        configfile : str
        """
 

        # job identifier
        self.name = "name"            
        # directory where configuration file is found.
        self.configdir = configdir    
        # configuration file for HYSPLIT runs
        self.configfile = configfile  
        # where to look for volcat data
        self.volcatdir = vdir        
        # volcano location
        self.vloc = [0, 0]         

        # list of subdirectories
        self.subdirlist = []

        # instance of the InverseAshEns class. 
        # InverseAshEns class maintains a list of InverseAsh objects.
        # Each InverseAsh object corresponding to one run in the ensemble.
        # one run meaning one set of unit source runs for the inversion.

        self.invens = ai.InverseAshEns(
            tdirlist,
            fnamelist,
            vdir,
            vid,
            configdir=configdir,
            configfile=configfile,
            verbose=False,
        )

        # particle size. TODO - more sophisticated particle size handling.
        self.invens.add_phash({"p060": 1})
        # self.shash contains information from the configuration file 
        self.shash = {}
        self.set_sourcehash()
        # directory information is obtained from self.shash, the config file.
        self.set_directories()
        # sets volcatdir
        self.set_vdir(vdir)
        # get the volcano location 
        self.vloc = self.get_vloc()  

    def get_vloc(self):
        chash = self.shash
        if "latitude" in chash.keys():
            latitude = chash["latitude"]
        else:
            latitude = None
        if "longitude" in chash.keys():
            longitude = chash["longitude"]
        else:
            longitude = None
        vloc = [longitude, latitude]
        return vloc

    def change_config_info(self,configdir,configfile):
        self.invens.invlist[0].add_config_info(configdir,configfile)
        self.set_sourcehash()
        self.set_directories()

    def set_sourcehash(self):
        self.shash = self.invens.invlist[0].inp

    def set_directories(self):
        wdir = self.shash["WORK_DIR"]
        hysplitdir = self.shash["HYSPLIT_DIR"]
        datadir = self.shash["WORK_DIR"]
        # place where inverse executable is found.
        execdir = self.shash['INV_DIR']
        self.invens.set_directory(wdir, execdir, datadir, hysplitdir)

    def set_vdir(self, vdir):
        self.invens.vdir = vdir
        for inv in self.invens.invlist:
            inv.vdir = vdir

    def add_inputs(self, inp):
        return -1

    def prepare_data(self, verbose=False):
        """ """
        for drange in self.create_dlist():
            if verbose:
                print(drange)
            self.invens.prepare_one_time(drange, verbose=verbose)
        # self.invens.make_tcm_mult()

    def time_lagged_tcms(
        self,
        tag,
        tlist,
        remove_cols=True,
        remove_rows=False,
        remove_sources=False,
        remove_ncs=0,
    ):
        """
        tiilist : list  of list of integers
        """
        for tiilist in tlist:
            self.create_tcm(
                tag,
                tiilist,
                remove_cols=remove_cols,
                remove_rows=remove_rows,
                remove_sources=remove_sources,
                remove_ncs=remove_ncs,
            )

    def create_tcm(
        self,
        tag,
        tiilist,
        remove_cols=True,
        remove_rows=False,
        remove_sources=False,
        remove_ncs=0,
        overwrite=False,
    ):
        """
        tiilist : list of integers
        remove_cols : boolean : remove columns with no model values
        remove_rows : boolean : remove clear sky observations (rows)
        remove_sources :
        remove_ncs : integer : remove clear sky observations which are close to plume.
        """
        done=False
        runtag = ai.create_runtag(
            tag, tiilist, remove_cols, remove_rows, remove_sources, remove_ncs
        )
        print("TAG", runtag)
        self.invens.set_subdirectory(runtag)
        subdir = self.invens.subdir
        self.add_subdir(subdir)

        enslist = self.invens.taglist
        tname = "{}.txt".format(tag)
        # check that the Parameters in file is where it is supposed to be.
        if not os.path.isfile(os.path.join(self.invens.wdir,'Parameters_in.dat.original')):
           pname = 'Parameters_in.dat.original'
           print('Parameters in file does not exist {}'.format(subdir))
           original = os.path.join(self.shash['DATA_DIR'],pname)
           if not os.path.isfile(original):
              print('Cannot find file {}'.format(original))
              return False 
           new = os.path.join(self.invens.wdir,pname)
           Helper.move(original,new)

        # check to see if the TCM has been created by checking for the out.dat file.
        for ens in enslist: 
            qfile = os.path.join(subdir,'{}_out.dat'.format(ens))
            if os.path.isfile(qfile):
               done=True
            else:
               done=False
        if not done or overwrite:
            # this creates the InverseAsh class tcm attribute.
            self.invens.make_tcm_mult(
                tiilist,
                remove_cols=remove_cols,
                remove_rows=remove_rows,
                remove_sources=remove_sources,
                remove_ncs=remove_ncs,
            )
            tcmstr = self.invens.write_tcm(os.path.join(subdir, tname), verbose=True)
      
            self.invens.run_tcm()
            self.invens.save_emis(runtag + '.csv', ens=False)
        # still need to make_tcm to get the self.tcm_columns information.
        else:
            print('ZZZZZZZZZZZZZZZZZZZZZZZZZZ')
            self.invens.make_tcm_mult(
                tiilist,
                remove_cols=remove_cols,
                remove_rows=remove_rows,
                remove_sources=remove_sources,
                remove_ncs=remove_ncs,
            )

        # this creates the emit-times file, control file and setup file.
        # self.invens.make_efile(vloc=self.vloc)

    def add_subdir(self,subdir):
        if subdir not in self.subdirlist:
           self.subdirlist.append(subdir)

    def plot_outdat(self):
        for subdir in self.subdirlist:
            print(subdir)
            self.invens.set_subdirectory(subdir.split('/')[-1])
            ilist = self.invens.read_outdat(eii=None)
            ilist[0].get_conc()
            try:
                ilist[0].plot_conc()
            except Exception as eee:
                print('Error in plotting {}'.format(eee))
                continue
            plt.show()


    def get_subdirlist(self,subdirlist):
        if isinstance(subdirlist,str):
           subdirlist=[subdirlist]
        if not isinstance(subdirlist,list):
           subdirlist = self.subdirlist
        return subdirlist

    def setup_hysplit_runs(self,subdirlist=None):

        subdirlist = self.get_subdirlist(subdirlist)

        # need to check that output from inversion exists
        def checkfiles(sdir):
            check = True
            return check
         
        for subdir in subdirlist:
            print('setting up in {}'.format(subdir))
            self.invens.subdir = subdir
            self.invens.make_efile(self.vloc)

    def run_hysplit(self, subdirlist=None):
        """
        list of subdirectories to run HYSPLIT in.
        ASCDATA.CFG, CONTROL, SETUP, EMITIMES files must already be created.
        """

        subdirlist = self.get_subdirlist(subdirlist)

        def checkfiles(sdir):
            check = True
            control_names = glob.glob(os.path.join(sdir, "CONTROL.*"))
            if len(control_names) < 1:
                check = False
            setup_names = glob.glob(os.path.join(sdir, "SETUP.*"))
            if len(control_names) < 1:
                check = False
            return check

        for subdir in subdirlist:
            if not checkfiles(subdir):
                print("CONTROL or SETUP file not found {}".format(subdir))
                continue
            cdump_names = glob.glob(os.path.join(subdir, "*cdump*"))
            if len(cdump_names) > 0:
                print("cdump files already found ", cdump_names)
                continue
            self.invens.subdir = subdir
            self.invens.run_hysplit()

    def create_cdump_netcdf(self, subdirlist,outname,verbose=False):
        subdirlist = self.get_subdirlist(subdirlist)
        netcdf_name = outname
        blist = []
        for iii, sdir in enumerate(subdirlist):
            source = sdir.split('/')
            source = source[-1]
            self.invens.subdir = sdir
            cdump_names = glob.glob(os.path.join(sdir, "*cdump*"))
            for cdump in cdump_names:
                ens = cdump.split('.')
                ens = ens[-1]
                source2 = source + '_' +  ens
                print('working on ', source)
                blist.append([cdump, "metfile", source2])
        print("Creating", blist)
        dset = hysplit.combine_dataset(blist, verbose=verbose)
        self.forecast = dset
        # creates a netcdf file.
        try:
            hysplit.write_with_compression(dset,outname)
        except Exception as eee:
            print('Could not produce netcdf file {}'.format(eee))
        #self.create_asheval(dset)
        return dset

    def create_asheval(self,model,verbose=False):
        tdir = self.shash['WORK_DIR']
        configdir = self.configdir
        configfile = self.configfile
        ensdim = 'ens'
        volcatdir = self.volcatdir
        vid = None
        self.aeval = ae.AshEval(tdir,model,volcatdir,vid,
                                configdir,configfile,verbose,ensdim)
    def run_asheval(self):
        for drange in self.create_dlist():
            self.aeval.prepare_one_time(drange)

    def create_dlist(self):
        """
        Uses the start date and duration of simulation from the
        inversion configuration file.
        """
        dlist = []
        sdate = self.shash["start_date"]
        dhr = self.shash["durationOfSimulation"]
        dt = datetime.timedelta(hours=1)
        for iii in np.arange(0, dhr):
            drange = [sdate, sdate + dt]
            dlist.append(drange)
            sdate += dt
        return dlist


def example():

    tdir = '/hysplit-users/alicec/tmp/Loa/RunA/'
    tdirlist = [tdir]
    volcatdir = '/pub/ECMWF/JPSS/VOLCAT/Files/Mauna_Loa/pc_corrected/'
    vid=None
    configdir = '/hysplit-users/alicec/utilhysplit/ashapp/'
    configfile = 'config.LoaA.txt'
    fnamelist = ['xrfile.invLoaA.nc']
    pai = RunInversion(tdirlist,fnamelist,volcatdir,vid,configdir,configfile)

    print(pai.invens.taglist)
    pai.invens.print_directories()

    # STEP 1 : run the inversion algorithm.     
    pai.prepare_data(verbose=True) 
    tii = 0
    pai.invens.compare_plotsA(daterange=None,tii=tii,zii=None)

    # STEP 2 : HYSPLIT runs
    #     

if __name__ == "__main__":
   setup_logger()
   if len(sys.argv) != 3:
      print_usage()
      sys.exit(1)


    



