import glob
import os

from utilvolc import qva_logic, volcat


class TestInputs:
    def __init__(self):
        inp = {}
        inp["JPSS_DIR"] = "/pub/jpsss_upload"
        inp["VOLCAT_LOGFILES"] = "/pub/ECMWF/JPSS/VOLCAT/LogFiles/"
        inp["VOLCAT_DIR"] = "/pub/ECMWF/JPSS/VOLCAT/Files/"
        self.inp = inp


def test_volcat_001():
    """
    test that VOLCAT netcdf files can be opened and certain data variables are present.
    """
    vname = "/Popocatepetl/"
    tinp = TestInputs()
    fdir = tinp.inp["VOLCAT_DIR"]
    fnames = glob.glob(fdir + vname + "*nc")
    # print(fdir + vname + "*nc")
    # print(fnames)
    # most recent file
    fname1 = max(fnames, key=os.path.getctime)
    # oldest file
    fname2 = min(fnames, key=os.path.getctime)
    for fname in [fname1, fname2]:
        dset = volcat.open_dataset(
            fname,
            gridspace=None,
            correct_parallax=None,
            mask_and_scale=True,
            decode_times=False,
        )
        assert "pc_latitude" in dset.data_vars
        assert "pc_longitude" in dset.data_vars
        assert "ash_mass_loading" in dset.data_vars
        assert "ash_cloud_height" in dset.data_vars
        assert "ash_mass_loading_total_mass" in dset.data_vars
        assert "feature_area" in dset.data_vars
        assert "effective_radius_of_ash" in dset.data_vars


def test_summaryfile_001():
    """
    test that summary files can be opened
    """
    tinp = TestInputs()
    fdir = tinp.inp["JPSS_DIR"]
    fnames = glob.glob(fdir + "/*json")
    # most recent file
    fname1 = max(fnames, key=os.path.getctime)
    # oldest file
    fname2 = min(fnames, key=os.path.getctime)
    for fname in [fname1, fname2]:
        sumf = qva_logic.SummaryFile(fname, fdir)
        datadf = sumf.open_dataframe()
        assert 1 == 1


def test_summaryfile_002():
    """
    test that summary files can be opened and contains correct columns.
    """
    tinp = TestInputs()
    fdir = tinp.inp["JPSS_DIR"]
    sumdf = qva_logic.get_summary_file_df(fdir)
    columns = sumdf.columns
    cnames = ["volcano_name", "volcano_gvp_id", "volcano_lat", "volcano_lon"]
    cnames.append("vaac_region")
    cnames.append("event_type")
    cnames.append("log_url")
    cnames.append("log")
    cnames.append("event_type")
    cnames.append("event date")
    cnames.append("event vid")
    cnames.append("event gid")
    cnames.append("summary date")
    cnames.append("summary file")
    cnames.append("satellite")
    cnames.append("volcano_elevation")
    cnames.append("volcano_country")
    for cnn in cnames:
        print(cnn, cnn in columns)
        assert cnn in columns


def test_eventfile_001():
    """
    test that event files can be opened and data has correct columns.
    """
    tinp = TestInputs()
    fdir = tinp.inp["VOLCAT_LOGFILES"]
    fnames = glob.glob(fdir + "/*json")
    # most recent file
    fname1 = max(fnames, key=os.path.getctime)
    # oldest file
    fname2 = min(fnames, key=os.path.getctime)
    for fname in [fname1, fname2]:
        efile = qva_logic.EventFile(fname, fdir)
        assert efile.open()
        columns = efile.df.columns
        cnames = ["sensor_name", "observation_date", "sensor_mode", "sensor_wmo_id"]
        cnames.append("start_coverage_time")
        cnames.append("observation_time")
        cnames.append("observation_epoch")
        cnames.append("feature_id")
        cnames.append("event_file")
        cnames.append("event_url")
        for cnn in cnames:
            print(cnn, cnn in columns)
            assert cnn in columns


test_volcat_001()
