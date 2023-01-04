import pandas as pd

def fixlondf(df, colname="lon", neg=True):
    """
    df : pandas DataFrame
    colname : str. name of column with longitude values to be converted.
    neg : boolean

    if neg=False
    convert longitude to degrees east (all positive values)

    if neg==True
    convert longitude to degrees west (all negative values).

    """
    if not neg:
        df[colname] = df.apply(lambda row: fixlon(row[colname]), axis=1)
    else:
        df[colname] = df.apply(lambda row: fixlon(row[colname]) - 360, axis=1)
    return df


def fixlon(xvar):
    if xvar < 0:
        return 360 + xvar
    else:
        return xvar


