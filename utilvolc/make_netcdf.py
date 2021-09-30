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
        d0: datetime object (date range from d1)
        sample_time_stamp: determining if time stamp is attributed to the end or beginning
        mettag: met field data tag (default: 'gfs0p25')
        volcid: volcano ID (string)
        volcname: volcano name (string)
        -----------------
        Functions:
        DI_combine: creates xarray dataset of ensemble Data Insertion hysplit simulations
        Cyl_combine: creates xarray dataset of cylinder source hysplit simulations
        Line_combine: creates xarray dataset of line source hysplit simulations
        create_ens: creates ensemble netcdf file of DI, Cyl, Line datasets
        create_volcat: creates regridded volcat netcdf - same as hysplit grid
        calc_stats: creates statistics netcdf containint Brier Score, Pearson Correlations for each ensemble member
        """
        self.d1 = d1
        self.d0 = d0
        self.sample_time_stamp = sample_time_stamp
        self.mettag = mettag
        self.volcid = volcid
        self.volcname = volcname

    def DI_combine(self, all_files, start=None, end=None, tdelta=10, verbose=False):
        """ Creates xarray dataset of Data Insertion hysplit simulations
        all_files: list of all cdump files from Data Insertion runs
        start: datetime object(first time to be included in ensemble (default: d1 - 24 hours))
        end: datetime object (last time to be included in ensemble (default: d1 - tdelta))
        tdelta: time interval of files to include in minutes (default: 10)
        verbose: (boolean)
        """
        if start == None:
            start = self.d1 - timedelta(hours=24)
        if end == None:
            end = self.d1
        files = []
        sourcetag = []
        metdatatag = []

        # Adding Data Insertion cdump files
        while start < end:
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
            blist, drange=[self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=verbose)
        hxrd.attrs['Volcano Name'] = self.volcname
        hxrd.attrs['Volcano ID'] = self.volcid
        del(hxrd.attrs['Coordinate time description'])
        return hxrd

    def Cyl_combine(self, allfiles, match=''):
        """ Creates xarray dataset of Cylinder source hysplit simulation
        allfiles: cdump file name(s) - use glob to create array
        match: string for matching in filename
        """
        if len(allfiles) > 1:
            files = []
            sourcetag = []
            metdatatag = []
            for f in allfiles:
                if match in f:
                    files.append(f)
                    s = f.find('cdump.')+len('cdump.')
                    e = f.find('hrs_')
                    typec = f[s:e]
                    sourcetag.append('Cyl_'+typec+'hrs_emis')
                    metdatatag.append(self.mettag)
            blist = list(zip(files, sourcetag, metdatatag))
            hxrcyl = hysplit.combine_dataset(
                blist, drange=[self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)
        else:
            hxrcyl = hysplit.combine_dataset([(allfiles[0], 'CylinderSource', self.mettag)], drange=[
                self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)
            del(hxrcyl.attrs['Coordinate time description'])
        return hxrcyl

    def Line_combine(self, allfiles, match=''):
        """Creates xarray dataset of Line source hysplit simulation
        filename: line file name(s) - use glob to make array
        match: string for matching in filename
        """
        if len(allfiles) > 1:
            files = []
            sourcetag = []
            metdatatag = []
            for f in allfiles:
                if match in f:
                    files.append(f)
                    s = f.find('cdump.')+len('cdump.')
                    e = f.find('hrs_')
                    typel = f[s:e]
                    sourcetag.append('Line_'+typel+'hrs_emis')
                    metdatatag.append(self.mettag)
            blist = list(zip(files, sourcetag, metdatatag))
            hxrline = hysplit.combine_dataset(
                blist, drange=[self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)
        else:
            hxrline = hysplit.combine_dataset([(allfiles[0], 'LineSource', self.mettag)], drange=[
                                              self.d0, self.d1], century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)

        hxrl = hysplit.open_dataset(allfiles[0], drange=[self.d0, self.d1],
                                    century=2000, sample_time_stamp=self.sample_time_stamp, verbose=False)

        unitmass, mass63 = hysp_func.calc_MER(hxrl)

        hxrline.attrs['Volcano Latitude'] = hxrl.attrs['Starting Latitudes'][0]
        hxrline.attrs['Volcano Longitude'] = hxrl.attrs['Starting Longitudes'][0]
        hxrline.attrs['Volcano Vent (m)'] = hxrl.attrs['Starting Heights'][0]
        hxrline.attrs['Plume Height (m)'] = hxrl.attrs['Starting Heights'][-1]
        hxrline.attrs['Line and Cylinder Source Date'] = hxrl.attrs['Source Date'][0][0:8]
        hxrline.attrs['Line and Cylinder Source Time'] = hxrl.attrs['Source Date'][0][9:16]
        hxrline.attrs['Mass Eruption Rate - Mastin'] = unitmass
        hxrline.attrs['Fine Ash MER - Mastin'] = mass63
        hxrline.attrs['MER Units'] = 'g/hr'
        del(hxrline.attrs['Coordinate time description'])
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
            return print(netdir+newfile+' created!')
        return hxr

    def create_volcat(self, volcatdir, netdir, volcdir, ensname, skipna=False, convert_nans=False, write=False):
        """ Creates netcdf of VOLCAT, gridded to HYSPLIT grid
        Inputs:
        volcatdir: directory of VOLCAT data
        netdir: directory of ensemble netcdf
        volcdir: directory of regridded VOLCAT data
        ensname: ensemble netcdf name
        skipna: (boolean) skips nan values when calculating average (default=False)
        convert_nans: (bolean) converts nan values to 0. when calculating average (default=False)
        write: (boolean) if True, netcdf file is written (default = False)
        Outputs:
        dnew: volcat xarray 
        if write: returns location of netcdf file"""

        from glob import glob
        from utilvolc import volcat

        # Will adjust this to use the VolcatName class in volcat.py
        vnames = glob(volcatdir+'*'+self.d0.strftime('s%Y%j_%H')+'*'+self.volcid+'*')
        vname2 = glob(volcatdir+'*'+self.d1.strftime('s%Y%j_%H%M%S')+'*'+self.volcid+'*')
        if len(vname2) != 0:
            vnames.append(vname2[0])
        dset = []
        x = 0
        while x < len(vnames):
            dset.append(volcat.open_dataset(vnames[x], decode_times=True, mask_and_scale=True))
            x += 1

        hxr = xr.open_dataset(netdir+ensname)
        # Regridding volcat files to hysplit lat/lon
        # dnew is xarray (time, x, y) with a regridded file for each time period
        dnew = volcat.average_volcat_new(dset, hxr, skipna=skipna, convert_nans=convert_nans)
        if write:
            # Removing old netcdf file if it exists
            volcnc = 'regridded_volcat_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
            if os.path.exists(volcdir+volcnc):
                os.remove(volcdir+volcnc)
            # Creating netcdf file
            dnew.to_netcdf(volcdir+volcnc)
            return print(volcdir+volcnc+' created!')
        else:
            return dnew

    def create_volcat_old(self, volcatdir, netdir, volcdir, write=False):
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
        if len(vname2) != 0:
            vdset = volcat.open_dataset(vname2[0], decode_times=True)
            dset.append(vdset)
            print(vname2[0])

        ensfile = 'ensemble_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'
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
        mass_avg.attrs['units'] = dset[0].ash_mass_loading.attrs['units']
        mass_avg.attrs['long_name'] = 'Average total column loading of ash in the highest continuous ash layer for the previous hour'
        mass_avg.attrs['FillValue'] = 'nan'
        hgt_max2 = hgt_max * 1000.  # Converting from km to m
        hgt_max2.attrs['long_name'] = 'Maximum cloud top height of the highest continuous ash layer for the previous hour'
        hgt_max2.attrs['units'] = 'm'
        hgt_max2.attrs['FillValue'] = 'nan'

        # Regridding volcat to hysplit resolution
        if len(vname2) != 0:
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

    def calc_stats(self, ensdir, volcdir, statdir, threshold, deltaz=1000., namestr='', pstr='par006', write=False):
        """ Calculates Brier Scores and Pearson Correlations for each ensemble
        member. Must convert concentration to mass loading using deltaz
        value and summing along z axis for comparison to VOLCAT data.

        IN PROGRESS - NEED TO MODIFY FOR USE WITH ENS DIMENSION 
        AND MULTIPLE VARIABLES
        Inputs:
        ensdir: ensemble netcdf directory
        volcdir: regridded volcat netcdf directory
        statdir: statistics netcdf directory
        threshold: list of floats (thresholds to calculate scores for)
        dimension: (string) dimension name to calculate statistics (source or ens)
        deltaz: z-axis interval size in meters (float) default: 1000.0 meters
        namestr: string in name for particular file
        pstr: string particle name (variable) 
        write: (boolean) whether or not to write netcdf file
        Outputs:
        If write = True, returns stats netcdf file
        if write = False, xarray dataset of statistical variables
        """
        from utilhysplit.evaluation import plume_stat as ps
        import numpy as np

        ensfile = 'ensemble_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+namestr+'_'+pstr+'.nc'
        volcfile = 'regridded_volcat_'+self.volcname+'_'+self.d1.strftime('%Y%m%d.%H%M%S')+'.nc'

        hxr = xr.open_dataset(ensdir+ensfile)
        vxr = xr.open_dataset(volcdir+volcfile)
        attrs = hxr.attrs
        # Summing along z makes hxr2 have units of g/m^2
        var = 'p'+pstr[3:]
        hxr2 = hxr[var] * deltaz
        hxr2 = hxr2.sum(dim='z')
        # Calculating number of forecasts with data in each grid box
        numfiles = len(hxr2.source)
        # Creating dummy xarray for merging below
        data = np.zeros((numfiles))
        statsxr = xr.DataArray(name='dummy', data=data,
                               dims='source', coords=[hxr2.source.values])
        # Calculations for various thresholds
        threshattr = []
        t = 0
        while t < len(threshold):
            # Converting to VOLCAT binary field for BS calculation
            if threshold[t] == 0.:
                ashmass = xr.where(vxr.ash_mass_loading.isel(time=-1) > threshold[t], 1., 0.)
                ashmassavg = xr.where(vxr.ash_mass_avg.isel(time=-1) >= threshold[t], 1., 0.)
            else:
                ashmass = xr.where(vxr.ash_mass_loading.isel(time=-1) >= threshold[t], 1., 0.)
                ashmassavg = xr.where(vxr.ash_mass_avg.isel(time=-1) >= threshold[t], 1., 0.)
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
                stats = ps.CalcScores(vxr.ash_mass_loading.isel(time=-1),
                                      hxr2.isel(source=a), threshold=threshold[t], verbose=False)
                PCcent, PCuncent = stats.calc_pcorr()
                PClistcent.append(PCcent.values)
                PClistuncent.append(PCuncent.values)

                stats2 = ps.CalcScores(vxr.ash_mass_avg.isel(
                    time=-1), hxr2.isel(source=a), threshold=threshold[t], verbose=False)
                PCcentavg, PCuncentavg = stats2.calc_pcorr()
                PClistcentavg.append(PCcentavg.values)
                PClistuncentavg.append(PCuncentavg.values)

                # Creating binary field of hysplit output
                hxr3 = xr.where(hxr2.isel(source=a) >= threshold[t], 1., 0.)
               # Calculating the BS values of each ensemble member
                BS = ps.calc_bs(hxr3, ashmass)
                BSlist.append(BS.values)
                BSavg = ps.calc_bs(hxr3, ashmassavg)
                BSlistavg.append(BSavg.values)
                a += 1
            # Adding Brier Scores to the netcdf, with dimension source
            thresh = str(threshold[t])
            threshattr.append(thresh)
            BSxr = xr.DataArray(BSlist, dims='source').load().rename('BS'+thresh)
            BSxr.attrs['long name'] = 'Brier Score compared to volcat'
            BSxr.attrs['threshold'] = thresh+' g/m^s'
            BSavgxr = xr.DataArray(BSlistavg, dims='source').load().rename('BSavg'+thresh)
            BSavgxr.attrs['long name'] = 'Brier Score compared to 1hr avg volcat'
            BSavgxr.attrs['threshold'] = thresh+' g/m^s'
            PCxr = xr.DataArray(PClistcent, dims='source').load().rename('PC'+thresh)
            PCxr.attrs['long name'] = 'Pattern Correlation (centered) compared to volcat'
            PCxr.attrs['threshold'] = thresh+' g/m^s'
            PCxruc = xr.DataArray(PClistuncent, dims='source').load().rename('PCuc'+thresh)
            PCxruc.attrs['long name'] = 'Pattern Correlation (uncentered) compared to volcat'
            PCxruc.attrs['threshold'] = thresh+' g/m^s'
            PCavgxr = xr.DataArray(PClistcent, dims='source').load().rename('PCavg'+thresh)
            PCavgxr.attrs['long name'] = 'Pattern Correlation (centered) compared to 1hr avg volcat'
            PCavgxr.attrs['threshold'] = thresh+' g/m^s'
            PCavgxruc = xr.DataArray(PClistuncent, dims='source').load().rename('PCucavg'+thresh)
            PCavgxruc.attrs['long name'] = 'Pattern Correlation (uncentered) compared to 1hr avg volcat'
            PCavgxruc.attrs['threshold'] = thresh+' g/m^s'
            statsxr = xr.merge([statsxr, BSxr, BSavgxr, PCxr, PCxruc, PCavgxr,
                                PCavgxruc], combine_attrs='drop_conflicts')
            t += 1
        # Dropping dummy variable
        statsxr = statsxr.drop(labels='dummy')
        statsxr.attrs['threshold'] = threshattr
        statsxr.attrs['threshold units'] = 'g/m^2'

        if write:
            # Removing and rewriting stats netcdf file
            statfile = self.volcname+'_statistics_'+self.d1.strftime('%Y%m%d.%H%M%S')+namestr+'_'+pstr+'.nc'
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
    from datetime import datetime
    from monetio.models import hysplit

    if MER == None:
        MER = hxrl.attrs['Fine Ash MER - Mastin']

    hxrcyl = hxrc * MER
    hxrline = hxrl * MER

    hxrcyl.attrs = hxrc.attrs
    hxrline.attrs = hxrl.attrs

    if MER != hxrl.attrs['Fine Ash MER - Mastin']:
        hxrline.attrs['User designated Fine Ash MER (Used)'] = MER

    hxr = xr.merge([hxrd, hxrcyl, hxrline], combine_attrs='drop_conflicts')
    # Need to confirm no nans in lat/lon arrays
    hxrnew = hysplit.reset_latlon_coords(hxr)
    return hxrnew


def calc_area(filename):
    from math import pi
    """Calculates the area (km^2) of each hysplit grid cell
    Converts degress to meters using a radius of 6378.137km.
    Input:
    filename: full file name of hysplit netcdf file
    output:
    area: xarray containing gridded area values
    """
    d2r = pi / 180.0  # convert degress to radians
    d2km = 6378.137 * d2r  # convert degree latitude to kilometers

    hxr = xr.open_dataset(filename)

    # Pulls out latitude and longitude
    lat = hxr.latitude
    lon = hxr.longitude
    latrad = lat * d2r  # Creating latitude array in radians
    coslat = np.cos(latrad) * d2km * d2km
    shape = np.shape(hxr.latitude)

    # Make shifted lat and shifted lon arrays to use for calculations
    lat_shift = lat[1:, :].values
    lon_shift = lon[:, 1:].values
    # Adding row/column of nans to shifted lat/lon arrays
    to_add_lon = np.empty([shape[0]]) * np.nan
    to_add_lat = np.empty([shape[1]]) * np.nan
    # Back to xarray for calculations
    lat2 = xr.DataArray(np.vstack((lat_shift, to_add_lat)), dims=['y', 'x'])
    lon2 = xr.DataArray(np.column_stack((lon_shift, to_add_lon)), dims=['y', 'x'])

    # area calculation
    area = abs(lat-lat2) * abs(abs(lon)-abs(lon2)) * coslat
    area.name = 'area'
    area.attrs['long_name'] = 'area of each hysplit lat/lon grid box'
    area.attrs['units'] = 'km^2'

    return area
