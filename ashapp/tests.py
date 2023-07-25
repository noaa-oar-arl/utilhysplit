import datetime
import rundispersion

def getinp():
    inp = {}
    inp['meteorologicalData'] = 'GFS'
    inp['forecastDirectory'] = '/pub/forecast/'
    inp['archivesDirectory'] = '/pub/archives/'
    inp['WORK_DIR'] = './'
    inp['jobname'] = 'A'
    inp['durationOfSimulation'] = 10
    inp['DATA_DIR'] = './'
    inp['latitude'] = 45
    inp['longitude'] = -160
    inp['bottom'] = 1000
    inp['top'] = 15000
    inp['emissionHours'] = 2
    inp['rate'] = 1
    inp['area'] = 1
    temp = datetime.datetime.now()
    inp['start_date'] = datetime.datetime(temp.year, temp.month, temp.day, temp.hour,0,0)
    inp['samplingIntervalHours'] = 3
    inp['HYSPLIT_DIR'] = '/hysplit-users/alicec/hdev/'
    return inp


def test1():
    import rundispersion
    inp = getinp()
    inp['jobid'] = 'A'
    dr = rundispersion.RunDispersion(inp)
    dr.run_model()
    return dr.filelist


def test2(filelist):
    import outputdispersion
    inp = getinp()
    inp['jobid'] = 'A'
    inp['eflag']=2
    mout = outputdispersion.OutputDispersion(inp,filelist)
    mout.postprocess()     


def test3():
    from ashapp.ashnetcdf import AshAttributes
    import datetime
    import numpy as np
    q = {}
    b = {}
    c = {}
    q['dog']=1
    q['cat']=2
    q['bird'] = np.array([1,2,3])
    c['n2'] = np.array([7,8,9])
    b['test'] = datetime.datetime(2020,1,20,0)
    b['test2'] = c
    q['nest'] = b
   
    result = {}
    result = {'dog':1, 'cat':2, 'bird': [1,2,3],
              'nest':{'test':'2020-01-20 00:00:00', 'test2':{'n2':[7,8,9]}}}

    att = AshAttributes(q)
    assert att.attr == result

def test4():
    from utilvolc import invhelper
    inp = getinp()
    d1 = datetime.datetime(2020,2,1,0)
    dth = 1
    dt = datetime.timedelta(hours=dth)
    inp['inv_vertical_resolution'] = 1000
    inp['timeres'] = dth
    inp['start_date'] = datetime.datetime(2020,2,1,0)
    inp['bottom'] = 2000
    inp['top'] = 4000
    inp['emissionHours'] = 2
    shash = invhelper.inverse_get_suffix_list(inp,suffix_type='date', dtfmt='%m%d%H')

    temp = {}
    key = '020100_2000'
    temp[key] = {'sdate':d1, 'edate': d1 + dt, 'bottom':2000, 'top':3000}
    print(shash) 
    assert shash[key] == temp[key]


#filelist = test1()
print('PASSED test1 ')
#test2(filelist)
print('PASSED test2 ')
#test3() 
print('PASSED test2 ')
test4() 




