import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt


def draw_map(fignum, fs=20):
    proj = ccrs.PlateCarree()
    fig = plt.figure(fignum, figsize=(16,9))
        
    ax = fig.add_subplot(1,1,1)
    ax = plt.axes(projection=proj)
    gl = ax.gridlines(draw_labels=True, linewidth = 0.2, color="gray")
    gl.labels_right=False
    gl.labels_top=False
    #ocean = cfeature.NaturalEarthFeature(category='physical',name='ocean',
    #        scale='50m', facecolor = cfeature.COLORE['water'])

    #ax.add_feature(states.edgecolor='gray')
    ax.add_feature(cfeature.LAND)
    return ax

