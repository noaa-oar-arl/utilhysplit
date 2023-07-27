
def fix_volc_name(volcname):
    """Fixes the volcano name if a comma, or space appear in the name"""
    if "," in volcname:
        s = volcname.find(",")
        tmp = volcname[:s]
        tmp2 = volcname[s + 2 :]
        volcname = tmp2 + "_" + tmp
    if " " in volcname:
        volcname = volcname.replace(" ", "_")
    return volcname

def get_latlon(data_dir):
    """Read csv file containing volcano name, latitude, longitude
    Inputs:
    data_dir: directory where Volcanoes.csv is located
    Outputs:
    volcdf: pandas dataframe of volcanos, latitude, longitude
    """
    import pandas as pd
    volcdf = pd.read_csv(data_dir + "Volcanoes.csv", sep=",")
    return volcdf
