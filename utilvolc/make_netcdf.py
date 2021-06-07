# make_netcdf.py
# Makes ensemble hysplit netcdf, volcat regrid netcdf, and stats netcdf
from datetime import datetime, timedelta
from monetio.models import hysplit
from utilvolc import hysp_func
import xarray as xr
import numpy as np
import os


class MakeNetcdf:
    def __init__(self, d1, d0, sample_time_stamp='end', mettag='gfs0p25', volcid='VolcanoID', volcname='VolcanoName'):
        """
        Class of tools to make netcdf files necessary for evaluation
        ----------------
        Inputs:
        d1: datetime object (time stamp of of ensemble)
        d0: datetime object (date range from d1 (default: d1 - 1 hour))
        sample_time_stamp: determining if time stamp is attributed to the end or beginning
        mettag: met field data tag (default: 'gfs0p25')
        volcid: volcano ID (string)
        volcname: volcano name (string)
        -----------------
        Functions:
        DI_combine: creates xarray dataset of ensemble Data Insertion hysplit simulations
        Cyl_combine: creates xarray dataset of cylinder source hysplit simulations
        Line_combine: creates xarray dataset of line source hysplit simulations
        Regrid_volcat:
        StatsBS: creates xarray dataarray of Brier Scores for each ensemble member
        StatsPC: creates xarray dataarray of Pearson Correlations for each ensemble member
        """
        self.d1 = d1
        self.d0 = d0
        self.sample_time_stamp = sample_time_stamp
        self.mettag = mettag
        self.volcid = volcid
        self.volcname = volcname

    def DI_combine(self, all_files, start=None, end=None, tdelta=10):
        """ Creates xarray dataset of Data Insertion hysplit simulations
        all_files: list of all cdump files from Data Insertion runs
        start: datetime object(first time to be included in ensemble (default: d1 - 24 hours))
        end: datetime object (last time to be included in ensemble (default: d1 - 2 hours))
        tdelta: time interval of files to include in minutes (default: 10)
        """
        if start == None:
            start = self.d1 - timedelta(hours=24)
        if end == None:
            end = self.d1 - timedelta(hours=2)
        files = []
        sourcetag = []
        metdatatag = []

        # Adding Data Insertion cdump files
        while start <= end:
            t = start.strftime('%Y%j_%H%M%S')
            for match in all_files:
                if t in match:
                    files.append(match)
                    sourcetag.append('DI_'+start.strftime('%Y%m%d.%H%M%S'))
                    metdatatag.append(self.mettag)
            start += timedelta(minutes=tdelta)
        # Creating tuple for hysplit.combine_dataset
        blist = list(zip(files, sourcetag, metdatatag))

        # Date Range for ensemble hour
        hxrd = hysplit.combine_dataset(
            blist, drange=[self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)
        hxrd.attrs['Volcano Name'] = self.volcname
        hxrd.attrs['Volcano ID'] = self.volcid
        return hxrd

    def Cyl_combine(self, directory, filename):
        """ Creates xarray dataset of Cylinder source hysplit simulation
        directory: directory of cylinder cdump file
        filename: cdump file name
        """
        hxrcyl = hysplit.combine_dataset([(directory+filename, 'CylinderSource', self.mettag)], drange=[
                                         self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)

        hxrc = hysplit.open_dataset(
            directory+filename, drange=[self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)

        hxrcyl.attrs['Cylinder Num Start Locations'] = hxrc.attrs['Number Start Locations']
        return hxrcyl

    def Line_combine(self, directory, filename):
        """Creates xarray dataset of Line source hysplit simulation
        directory: directory of line cdump file
        filename: line file name
        """
        hxrline = hysplit.combine_dataset([(directory+filename, 'LineSource', self.mettag)], drange=[
                                          self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)

        hxrl = hysplit.open_dataset(
            directory+filename, drange=[self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)

        unitmass, mass63 = hysp_func.calc_MER(hxrl)

        hxrline.attrs['Volcano Latitude'] = hxrl.attrs['Starting Locations'][0][0]
        hxrline.attrs['Volcano Longitude'] = hxrl.attrs['Starting Locations'][0][1]
        hxrline.attrs['Volcano Vent (m)'] = hxrl.attrs['Starting Locations'][0][2]
        hxrline.attrs['Plume Height (m)'] = hxrl.attrs['Starting Locations'][1][2]
        hxrline.attrs['Line and Cylinder Source Date'] = hxrl.attrs['Source Date'][0].strftime('%Y%m%d')
        hxrline.attrs['Line and Cylinder Source Time'] = hxrl.attrs['Source Date'][0].strftime('%H%M%S')
        hxrline.attrs['Line Num Start Locations'] = hxrl.attrs['Number Start Locations']
        hxrline.attrs['Mass Eruption Rate - Mastin'] = unitmass
        hxrline.attrs['Fine Ash MER - Mastin'] = mass63
        hxrline.attrs['MER Units'] = 'g/hr'
        return hxrline

    def create_ens(self, hxrd, hxrc, hxrl, netdir, MER=None, write=False):
        """ Creates netcdf file of merged dataarrays
        Inputs:
        hxrd: data array of DI hysplit simulations
        hxrc: data array of cylinder source hysplit simulations
        hxrl: data array of line source hysplit simulations
        netdir: directory of ensemble netcdf file (string)
        MER: float - default is Mastin MER from hxrl attributes
        write: (boolean) If true, netcdf file is written (default=False)
        Output:
        If write=True, merged netcdf file
        if write=False, xarray dataset of merged dataarrays"""

        hxr = ens_merge(hxrd, hxrc, hxrl, MER=MER)
        if write:
            newfile = 'ensemble_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
            if os.path.exists(netdir+newfile):
                os.remove(netdir+newfile)
                hxr.to_netcdf(netdir+newfile)
            return print(newfile+' created!')
        return hxr

    def create_volcat(self, volcatdir, netdir, volcdir, write=False):
        """ Creates netcdf of VOLCAT, gridded to HYSPLIT grid
        Inputs:
        volcatdir: directory of VOLCAT data
        netdir: directory of ensemble netcdf
        volcdir: directory of regridded VOLCAT data
        write: (boolean) if True, netcdf file is written (default = False)
        Outputs:
        If write = True, netcdf file of regridded VOLCAT
        if write = False, xarray dataset of regridded VOLCAT
        """
        from glob import glob
        from utilvolc import volcat

        vnames = glob(volcatdir+'*'+self.d0.strftime('%Y%j_%H')+'*'+self.volcid+'*')
        vname2 = glob(volcatdir+'*'+self.d1.strftime('%Y%j_%H%M%S')+'*'+self.volcid+'*')
        dset = []
        x = 0
        while x < len(vnames):
            dset.append(volcat.open_dataset(vnames[x], decode_times=True))
            x += 1
        if os.path.exists(vname2[0]):
            vdset = volcat.open_dataset(vname2[0], decode_times=True)
            dset.append(vdset)
        print(vname2[0])
        ensfile = 'ensemble_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
        # ensfile = 'ensemble_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
        hxr = xr.open_dataset(netdir+ensfile)

        avgmass, maxhgt = volcat.average_volcat(dset, hxr)
        mass_avg = avgmass.load().rename('volcat_AshMass_avg')
        hgt_max = maxhgt.load().rename('volcat_AshHgt_max')
        # Need to fix longitude - HYSPLIT are from 0-360., VOLCAT are from -180. to 180.
        if np.min(hgt_max.longitude) < 0.:
            longit = hgt_max.longitude
            longit = longit.where(longit > 0, longit + 360.)
            hgt_max['longitude'] = longit
            mass_avg['longitude'] = longit

        # Assigning attributes
        mass_avg.attrs['units'] = vdset.ash_mass_loading.attrs['units']
        mass_avg.attrs['long_name'] = 'Average total column loading of ash in the highest continuous ash layer for the previous hour'
        mass_avg.attrs['FillValue'] = 'nan'
        hgt_max2 = hgt_max * 1000.  # Converting from km to m
        hgt_max2.attrs['long_name'] = 'Maximum cloud top height of the highest continuous ash layer for the previous hour'
        hgt_max2.attrs['units'] = 'm'
        hgt_max2.attrs['FillValue'] = 'nan'

        # Regridding volcat to hysplit resolution
        if os.path.exists(vname2[0]):
            mass_now = hxr.monet.remap_nearest(volcat.get_mass(vdset))
            hgt_now = hxr.monet.remap_nearest(volcat.get_height(vdset))
            # Renaming data variables and merging into dataset
            mass_now = mass_now.load().rename('volcat_AshMass')
            hgt_now = hgt_now.load().rename('volcat_AshHgt')
            # Need to fix longitude - HYSPLIT are from 0-360., VOLCAT are from -180. to 180.
            if np.min(hgt_now.longitude) < 0.:
                longit = hgt_now.longitude
                longit = longit.where(longit > 0, longit + 360.)
                hgt_now['longitude'] = longit
                mass_now['longitude'] = longit

            # Assigning attributes
            mass_now.attrs['long_name'] = vdset.ash_mass_loading.attrs['long_name']
            mass_now.attrs['units'] = vdset.ash_mass_loading.attrs['units']
            mass_now.attrs['FillValue'] = 'nan'
            hgt_now2 = hgt_now * 1000.  # Converting from km to m
            hgt_now2.attrs['long_name'] = vdset.ash_cloud_height.attrs['long_name']
            hgt_now2.attrs['FillValue'] = 'nan'
            hgt_now2.attrs['units'] = 'm'

            # Assigning coordinates for avg/max arrays
            mass_avg = mass_avg.expand_dims(dim='time').assign_coords(time=mass_now.time)
            hgt_max2 = hgt_max2.expand_dims(dim='time').assign_coords(time=hgt_now2.time)
            # Merging datasets
            hxrnew = xr.merge([mass_now, hgt_now2, mass_avg, hgt_max2], combine_attrs='drop_conflicts')
        else:
            hxrnew = xr.merge([mass_avg, hgt_max2], combine_attrs='drop_conflicts')
        attrs = dset[-1].attrs
        hxrnew.attrs = attrs
        hxrnew.attrs['starting time for avg'] = str(dset[0].time.values[0])
        hxrnew.attrs['ending time for avg'] = str(dset[-1].time.values[0])
        hxrnew.attrs['volcano ID'] = self.volcid

        if write:
            # Removing old netcdf file if it exists
            volcnc = 'regridded_volcat_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
            if os.path.exists(volcdir+volcnc):
                os.remove(volcdir+volcnc)
            # Creating netcdf file
            hxrnew.to_netcdf(volcdir+volcnc)
            return print(volcnc+' created!')
        return hxrnew

    def copy_keys(self, xra1, xra2):
        """Copies specific attributes from VOLCAT xarray to another xarray
        Inputs:
        xra1: array to take attributes from
        xra2: array to copy attributes to
        Returns: xra2"""

        xra2.attrs['title'] = xra1.attrs['title']
        xra2.attrs['institution'] = xra1.attrs['institution']
        xra2.attrs['volcano_name'] = self.volcname
        xra2.attrs['volcano_ID'] = self.volcid
        xra2.attrs['instrument_ID'] = xra1.attrs['instrument_ID']
        xra2.attrs['from_file'] = xra1.attrs['dataset_name']
        return xra2

    def calc_stats(self, ensdir, volcdir, statdir, threshold, deltaz=1000., write=False):
        """ Calculates Brier Scores and Pearson Correlations for each ensemble
        member. Must convert concentration to mass loading using deltaz
        value and summing along z axis for comparison to VOLCAT data.
        Inputs:
        ensdir: ensemble netcdf directory
        volcdir: regridded volcat netcdf directory
        statdir: statistics netcdf directory
        threshold: list of floats (thresholds to calculate scores for)
        deltaz: z-axis interval size in meters (float) default: 1000.0 meters
        write: (boolean) whether or not to write netcdf file
        Outputs:
        If write = True, returns stats netcdf file
        if write = False, xarray dataset of statistical variables
        """
        from utilhysplit.evaluation import plume_stat as ps
        import numpy as np

        ensfile = 'ensemble_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
        #ensfile = 'ensemble_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
        volcfile = 'regridded_volcat_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
        #volcfile = 'regridded_volcat_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'

        hxr = xr.open_dataset(ensdir+ensfile).squeeze()
        vxr = xr.open_dataset(volcdir+volcfile).squeeze()
        attrs = hxr.attrs
        # Summing along z makes hxr2 have units of g/m^2
        hxr2 = hxr.p006 * deltaz
        hxr2 = hxr2.sum(dim='z').squeeze()
        # Calculating number of forecasts with data in each grid box
        numfiles = len(hxr2.source)
        # Creating dummy xarray for merging below
        data = np.zeros((numfiles))
        statsxr = xr.DataArray(name='dummy', data=data, attrs=attrs,
                               dims='source', coords=[hxr2.source.values])
        # Calculations for various thresholds
        t = 0
        while t < len(threshold):
            # Converting to VOLCAT binary field for BS calculation
            ashmass = xr.where(vxr.volcat_AshMass.squeeze() >= threshold[t], 1., 0.)
            ashmassavg = xr.where(vxr.volcat_AshMass_avg.squeeze() >= threshold[t], 1., 0.)
            # Calculating Brier Score of each ensemble member
            a = 0
            BSlist = []
            BSlistavg = []
            PClistcent = []
            PClistuncent = []
            PClistcentavg = []
            PClistuncentavg = []
            while a < numfiles:
                # Calculating pattern correlation coefficients, centered and uncentered
                stats = ps.CalcScores(vxr.volcat_AshMass.squeeze(),
                                      hxr2[a, :, :], threshold=threshold[t], verbose=False)
                PCcent, PCuncent = stats.calc_pcorr()
                PClistcent.append(PCcent.values)
                PClistuncent.append(PCuncent.values)

                stats2 = ps.CalcScores(vxr.volcat_AshMass_avg.squeeze(),
                                       hxr2[a, :, :], threshold=threshold[t], verbose=False)
                PCcentavg, PCuncentavg = stats2.calc_pcorr()
                PClistcentavg.append(PCcentavg.values)
                PClistuncentavg.append(PCuncentavg.values)

                # Creating binary field of hysplit output
                hxr3 = xr.where(hxr2[a, :, :] >= threshold[t], 1., 0.)
               # Calculating the BS values of each ensemble member
                BS = ps.calc_bs(hxr3, ashmass)
                BSlist.append(BS.values)
                BSavg = ps.calc_bs(hxr3, ashmassavg)
                BSlistavg.append(BSavg.values)
                a += 1
            # Adding Brier Scores to the netcdf, with dimension source
            if threshold[t] == 0.1:
                thresh = '0p1'
            else:
                thresh = str(threshold[t])

            threshstr = str(threshold[t])+' g/m^2'
            BSxr = xr.DataArray(BSlist, dims='source').load().rename('BS'+thresh)
            BSxr.attrs['long name'] = 'Brier Score compared to volcat'
            BSxr.attrs['threshold'] = threshstr
            BSavgxr = xr.DataArray(BSlistavg, dims='source').load().rename('BSavg'+thresh)
            BSavgxr.attrs['long name'] = 'Brier Score compared to 1hr avg volcat'
            BSavgxr.attrs['threshold'] = threshstr
            PCxr = xr.DataArray(PClistcent, dims='source').load().rename('PC'+thresh)
            PCxr.attrs['long name'] = 'Pattern Correlation (centered) compared to volcat'
            PCxr.attrs['threshold'] = threshstr
            PCxruc = xr.DataArray(PClistuncent, dims='source').load().rename('PCuc'+thresh)
            PCxruc.attrs['long name'] = 'Pattern Correlation (uncentered) compared to volcat'
            PCxruc.attrs['threshold'] = threshstr
            PCavgxr = xr.DataArray(PClistcent, dims='source').load().rename('PC'+thresh)
            PCavgxr.attrs['long name'] = 'Pattern Correlation (centered) compared to 1hr avg volcat'
            PCavgxr.attrs['threshold'] = threshstr
            PCavgxruc = xr.DataArray(PClistuncent, dims='source').load().rename('PCuc'+thresh)
            PCavgxruc.attrs['long name'] = 'Pattern Correlation (uncentered) compared to 1hr avg volcat'
            PCavgxruc.attrs['threshold'] = threshstr

            statsxr = xr.merge([statsxr, BSxr, BSavgxr, PCxr, PCxruc, PCavgxr,
                                PCavgxruc], combine_attrs='drop_conflicts')
            t += 1
        # Dropping dummy variable
        statsxr = statsxr.drop(labels='dummy')

        if write:
            # Removing and rewriting stats netcdf file
            statfile = self.volcname+'_statistics_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
            if os.path.exists(statdir+statfile):
                os.remove(statdir+statfile)
            # Creating netcdf file
            statsxr.to_netcdf(statdir+statfile)
            return print(statfile+' created!')
        return statsxr


def ens_merge(hxrd, hxrc, hxrl, MER=None):
    """Merges the DI, cylinder, and line source dataarrays, multiplies the
    cylinder and line source data arrays by the MER
    hxrd: data array of DI hysplit simulations
    hxrc: data array of cylinder source hysplit simulations
    hxrl: data array of line source hysplit simulations
    MER: float - default is Mastin MER from hxrl attributes"""

    if MER == None:
        MER = hxrl.attrs['Fine Ash MER - Mastin']

    hxrcyl = hxrc * MER
    hxrline = hxrl * MER

    hxrcyl.attrs = hxrc.attrs
    hxrline.attrs = hxrl.attrs

    if MER != hxrl.attrs['Fine Ash MER - Mastin']:
        hxrline.attrs['User designated Fine Ash MER (Used)'] = MER

    hxr = xr.merge([hxrd, hxrcyl, hxrline], combine_attrs='drop_conflicts')
    return hxr
