import glob
import os
import datetime
from utilvolc import runhelper



class TestInputs:
    def __init__(self):
        inp = {}
        inp["JPSS_DIR"] = "/pub/jpsss_upload"
        inp["VOLCAT_LOGFILES"] = "/pub/ECMWF/JPSS/VOLCAT/LogFiles/"
        inp["VOLCAT_DIR"] = "/pub/ECMWF/JPSS/VOLCAT/Files/"
        self.inp = inp


def test_filenamecomposer_001():
    """
    """
    wdir = './'
    jobid = '999'
    jobname = 'Popo'
    fnc = runhelper.JobFileNameComposer(wdir,jobid,jobname)
    stagelist=[0,1,'a']
    for stage in stagelist:
        print('stage', stage)
        print('control filename',  fnc.get_control_filename(stage))
        print('setup filename',  fnc.get_setup_filename(stage))
        print('-----------------------')
        print('cdump filename',  fnc.get_cdump_filename(stage))
        print('pardump filename', fnc.get_pardump_filename(stage))
        print('tdump filename',  fnc.get_tdump_filename(stage))

        print('suffix',  fnc.make_suffix(stage))
        print(fnc.get_concplot_filename(stage))
        print(fnc.get_massloading_filename(stage,frame=0,ptype='gif'))
        print(fnc.get_massloading_filename(stage,frame=0,ptype='gif'))
        #print(fnc.get_trajplot_filename(stage,frame=0,ptype='gif'))
        #print(fnc.get_trajplot_filename(stage,frame=1,ptype='gif'))
        print('*********************')

    # testing ensemble files
    stage = 'gec00_102119_9880'
    target = 'CONTROL.{}_{}_{}'.format(jobname,jobid,stage)
    print(fnc.get_control_filename(stage))
    assert fnc.get_control_filename(stage)==target

    target = 'xrfile.{}_{}.nc'.format(jobname,jobid)
    test = fnc.get_xrfile()
    print('xrfile', test, target)
    assert target == test

    stage=1
    target = 'awips2.{}_{}_{}.nc'.format(jobname,jobid,stage)
    test = fnc.get_awips_filename(stage=stage)
    print(test, target)
    assert target == test

test_filenamecomposer_001()
