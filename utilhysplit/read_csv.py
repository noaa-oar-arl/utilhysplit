#read_csv.py
#Reads in csv files
#For use with MONET
import pandas as pd

#Opens CSV file and reads it into a Pandas DataFrame
def open_file(fname,delimiter):
    #print(fname)
    readCSV = pd.read_csv(fname, delimiter = delimiter, encoding = 'ISO-8859-1')
    return readCSV 

#Returns CSV file headers (assuming one line of header info)
def find_headers(readCSV):
    headers = list(readCSV.columns.values)
    return headers

