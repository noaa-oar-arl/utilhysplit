import numpy as np
import pandas as pd


def growth1(xo, yo, vx, vy, t=np.arange(0, 2)):
    """
    Calculate the growth of a rectangular area over time.
    
    This function computes how a rectangle with initial dimensions (dx, dy)
    grows over time when its sides expand at constant velocities (vx, vy).
    a = (xo+vx*t) * (yo+vy*t)
    a = xoyo + (vx*yo + vy*xo)*t + vx*vy*t^2
    The function returns the area of the rectangle and its dimensions at each time point.

    Parameters
    ----------
    dx : float
        Initial width of the rectangle.
    dy : float
        Initial height of the rectangle.
    vx : float
        Expansion velocity in the x-direction.
    vy : float
        Expansion velocity in the y-direction.
    t : numpy.ndarray, optional
        Time points at which to calculate the dimensions, 
        default is np.arange(0, 2).
    
    Returns
    -------
    a1 : numpy.ndarray
        Area of the rectangle at each time point.
    dx1 : numpy.ndarray
        Width of the rectangle at each time point.
    dy1 : numpy.ndarray
        Height of the rectangle at each time point.
    """
    ao = xo*yo + t*0 
    x1 = xo + vx*t
    y1 = yo + vy*t
    a1 = x1*y1
    df = pd.DataFrame({
        'ao': ao,
        'time': t,
        'area': a1,
        'lx': x1,
        'ly': y1
    })
    return df


def post_eruption(df,t=1,vx=1,vy=1):
    df2 = df.copy()
    df2 = df2[df2['time'] == t]
    dlist = []
    for row in df2.itertuples():
        xo = row.lx
        yo = row.ly
        dlist.append(growth1(xo, yo, vx, vy, t=np.arange(0, 24)))
    df3 = pd.concat(dlist, ignore_index=True)
    return df3 

def case1_eruption():
    a = np.array([0.4,1,2,5,10,16,24,50,121,274,619])
    yo = np.array([0.3,0.4,0.6,1.0,1.4,1.8,2.2,3.2,4.9,7.4,11.1])
    xo = np.array([1.4,2.1,3.2,4.8,7.2,8.9,10.9,15.9,24.6,37.0,55.6])

    ws = 70 #m/s
    ws = ws / 1000 # km/s
    #ws = ws * 60  # km/min
    vx = ws * 3600 # km/h

    # eruption lasts for one hour.
    # vx is the wind speed as ash is advected downwind.
    # vy is only due to dispersion so much smaller. 
    vy = vx * 0.01 # 1% of vx
    t = np.arange(0, 2, 1)  # time in hours
    alist = []
    for i in range(len(a)):
        df = growth1(xo[i], yo[i], vx, vy, t)
        alist.append(df)
    alist = pd.concat(alist, ignore_index=True)
    return alist

def rsd(df):
    """
    Calculate the relative standard deviation of the area over time.
    
    This function computes the relative standard deviation (RSD) of the area
    of the ash cloud at each time point, which is a measure of the variability
    of the area relative to its mean.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing 'ao', 'time', and 'area' columns.
    
    Returns
    -------
    rsd : pandas.DataFrame
        DataFrame with 'ao', 'time', and 'rsd' columns.
    """
    rsd = df.groupby('time').agg({'area': lambda x: np.std(x,ddof=1) / np.mean(x)}).reset_index()
    rsd.columns = ['time', 'rsd']
    return rsd

def plotdf(df,logy=False):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))
    for a in df['ao'].unique():
        sub_df = df[df['ao'] == a]
        plt.plot(sub_df['time'], sub_df['area'], label=f'Initial Area: {a:.2f}')
    if logy:
        plt.yscale('log') 
    plt.xlabel('Time (hours)')
    plt.ylabel('Area (km²)')
    plt.title('Growth of Ash Cloud Over Time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_with_custom_ticks(df, axis='y', logy=False):
    """
    Create a plot with custom tick marks at specified intervals.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the data to plot.
    axis : str, optional
        Which axis to customize ticks for ('x', 'y', or 'both'), default is 'y'.
    logy : bool, optional
        Whether to use logarithmic scale for y-axis, default is False.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))
    
    # Plot the data
    for a in df['ao'].unique():
        sub_df = df[df['ao'] == a]
        plt.plot(sub_df['time'], sub_df['area'], label=f'Initial Area: {a:.2f}')
    
    # Set custom tick positions
    custom_ticks = [0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5]
    
    # Apply custom ticks to specified axis
    if axis == 'x' or axis == 'both':
        plt.xticks(custom_ticks)
        plt.xlim([min(custom_ticks), max(custom_ticks)])
    
    if axis == 'y' or axis == 'both':
        plt.yticks(custom_ticks)
        plt.ylim([min(custom_ticks), max(custom_ticks)])
    
    # Format the tick labels to be more readable
    if axis == 'y' or axis == 'both':
        ax = plt.gca()
        ax.get_yaxis().set_major_formatter(
            plt.FuncFormatter(lambda x, loc: f"{int(x/1000)}k" if x > 0 else "0")
        )
    
    if logy:
        plt.yscale('log')
    
    plt.xlabel('Time (hours)')
    plt.ylabel('Area (km²)')
    plt.title('Growth of Ash Cloud Over Time')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
