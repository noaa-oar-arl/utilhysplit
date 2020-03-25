from timezonefinder import TimezoneFinder
import pandas as pd
import datetime

##adds time zone information and time offset to csv file with info about cems sites.

tf = TimezoneFinder()
fname = '/n-home/alicec/MONET/monet/data/cem_facility_loc.csv'

df = pd.read_csv(fname)

print(df[0:10])
print(df.columns.values)
print(df['longitude'])
print(df['latitude'])
def tzfind(x):
    return tf.timezone_at(lng=x['longitude'], lat=x['latitude'])

df['timezone'] = df.apply(tzfind, axis=1)

def getdate(x):
    t1 =  pd.Timestamp(datetime.datetime(2010,2,1,0)).tz_localize(x['timezone'])
    t2 =  t1.tz_convert('utc')
    t1 = t1.tz_localize(None)
    t2 = t2.tz_localize(None)
    return (t2 - t1).seconds/3600.0

df['time_offset'] = df.apply(getdate, axis=1)


print(df[0:10])
print(df['timezone'].unique())
print(df['time_offset'].unique())

fname2 = 'cemsinfo.csv'
df.to_csv(fname2, header=True)
