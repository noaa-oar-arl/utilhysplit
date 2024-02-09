import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import kasatochi_example as ke
from utilhysplit import par2conc
import seaborn as sns


def make_legend(keylist, chash):
    for tag in keylist:
        # print(chash[tag])
        plt.plot(
            1,
            1,
            marker="o",
            markeredgecolor=chash[tag],
            markerfacecolor=chash[tag],
            linewidth=0,
            label=tag.strftime("%Y %b %d %H:%M UTC"),
        )
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    plt.close()
    figlegend = plt.figure()
    ax1 = figlegend.add_subplot(1, 1, 1)
    ax1.legend(handles, labels, loc="center")
    ax1.axis("off")
    plt.savefig("feature_legend.png")
    plt.show()


def makelist(pdumpT, method="bgm", nnn=10, poll=4, lev=[1000, 20000]):
    # Using the GMM
    Trunlist = []
    datelist = []
    iii = 0
    for edate in pdumpT.date.unique():
        tempdf = ke.process_pdump(pdumpT, edate, lev=lev, poll=poll)
        if iii == 0:
            tfit = par2conc.par2fit(tempdf, method=method, nnn=nnn)
        else:
            tfit = par2conc.par2fit(tempdf, method=method, pfit=tfitp, nnn=nnn)
        tfitp = tfit
        iii += 1
        Trunlist.append(tfit)
        edate = datetime.datetime.strptime(str(edate), "%Y-%m-%dT%H:%M:%S.000000000")
        datelist.append(edate)
    return Trunlist, datelist


def process_fig(ax,dim,fname,add_text=None,xlim=None,ylim=None):
    ax.set_aspect('auto')
    if dim=='lat':
       ax.set_ylabel('height (km)')
    else:
       ax.set_ylabel('latitude')
    ax.set_xlabel('longitude')
    if xlim:
       print('setting xlim', xlim)
       ax.set_xlim(xlim[0],xlim[1])
    if ylim:
       ax.set_ylim(ylim[0],ylim[1])
    if add_text:
       iii = add_text[0]
       jjj = add_text[1]
       label = make_label(iii,jjj)
       dstr = add_text[2]
       ax.text(0.05,1.00,label+ ' ' + dstr,transform=ax.transAxes,bbox=dict(facecolor='white'))
       #ax.text(0.05,1.0,dstr,transform=ax.transAxes,bbox=dict(facecolor='white'))
    plt.savefig('{}.png'.format(fname))


def compare2(ft1, ft2, iii, dim="ht", name1="1", name2="2",
             cmap="plasma",fname='fig', add_text=False,
             xlim=None, ylim=None):
    sns.set(font_scale=1)
    sns.set_style('whitegrid')
    gfit1 = ft1.runlist[iii].gfit
    gfit2 = ft2.runlist[iii].gfit
    xr1 = ft1.runlist[iii].xra
    xr2 = ft2.runlist[iii].xra
    label = "score"
    dstr = ft1.datelist[iii].strftime("%Y %m/%d %H:%M UTC")
    dstr2 = ft1.datelist[iii].strftime("%m%d%H")
    print(dstr)
    print(name1, name1)
    par2conc.scatter(xr1, gfit1, cmap=cmap, labels=label, dim=dim)
    process_fig(plt.gca(), dim, fname+name1+name1+dstr2+dim, [name1,name1,dstr],
                xlim,ylim)
    plt.show()
    print(name2, name2)
    par2conc.scatter(xr2, gfit2, cmap=cmap, labels=label, dim=dim)
    process_fig(plt.gca(), dim, fname+name2+name2+dstr2+dim,[name2,name2,dstr],
                xlim,ylim)
    plt.show()
    print(name1, name2)
    par2conc.scatter(xr1, gfit2, cmap=cmap, labels=label, dim=dim)
    process_fig(plt.gca(), dim, fname+name1+name2+dstr2+dim,[name2,name1,dstr],
                xlim,ylim)
    plt.show()
    print(name2, name1)
    par2conc.scatter(xr2, gfit1, cmap=cmap, labels=label, dim=dim)
    process_fig(plt.gca(), dim, fname+name2+name1+dstr2+dim,[name1,name2,dstr],
                xlim,ylim)
    plt.show()

def make_label(iii, jjj):
    #return "{} with fit to {}".format(name2, name1)
    substr = '{}{}'.format(iii,jjj)
    return "S$_{" + substr + "}$"


def find_n(df1, df2, method='gmm',
           nlist=np.arange(5,80,5), 
           name1="1", name2="2",verbose=False):
    sns.set()
    sns.set_style('whitegrid')
    # input two pandas dataframes.
    score = {}
    score["1_2"] = []
    score["2_1"] = []
    score["2_2"] = []
    score["1_1"] = []
    for nnn in nlist:
        if verbose:
           print('working on {}'.format(nnn))
        # create fits for both DataFrames
        mfit1 = par2conc.par2fit(df1, method=method, nnn=nnn)
        mfit2 = par2conc.par2fit(df2, method=method, nnn=nnn)
        # compuate scores.
        s12, s21, s11, s22 = compare_sub(mfit1, mfit2)
        score["1_2"].append(s12)
        score["2_1"].append(s21)
        score["1_1"].append(s11)
        score["2_2"].append(s22)
    ax = plot_score(score, nlist, name1, name2)
    plt.xlabel("number of Gaussians")
    plt.tight_layout()
    return score   
 

def compare_sub(fit1, fit2):
    gfit1 = fit1.gfit
    gfit2 = fit2.gfit
    xr1 = fit1.xra
    xr2 = fit2.xra
    s12 = gfit1.score(xr2)
    s21 = gfit2.score(xr1)
    s11 = gfit1.score(xr1)
    s22 = gfit2.score(xr2)
    return s12, s21, s11, s22


def compare(ft1, ft2, name1="1", name2="2"):
    """
    ft1 FeatureTracker
    ft2 FeatureTracker
    """
    plt.figure(figsize=(10, 5))
    sns.set()
    sns.set(font_scale=2)
    sns.set_style("whitegrid")
    #sns.set(font_scale=2)
    # for now assume same dates.
    datelist = ft1.datelist

    score = {}
    score["1_2"] = []
    score["2_1"] = []
    score["2_2"] = []
    score["1_1"] = []

    for iii in np.arange(0, len(datelist)):
        s12, s21, s11, s22 = compare_sub(ft1[iii], ft2[iii])

        #gfit1 = ft1.runlist[iii].gfit
        #gfit2 = ft2.runlist[iii].gfit
        #xr1 = ft1.runlist[iii].xra
        #xr2 = ft2.runlist[iii].xra
        # use points in 2 with fit from 1
        score["1_2"].append(s12)
        score["2_1"].append(s21)
        score["1_1"].append(s11)
        score["2_2"].append(s22)
    ax = plot_score(score, datelist, name1, name2)
    plt.xlabel("Day Hour in August 2008 (UTC)")
    dform = DateFormatter("%d %H")
    ax.xaxis.set_major_formatter(dform)
    return score

def plot_score(score, xvals, name1, name2):
    clr = {}
    clr["1_1"] = "-ko"
    clr["2_2"] = "--g^"
    clr["1_2"] = "--r."
    clr["2_1"] = "--b."
    label = {}
    label["1_2"] = make_label(name1, name2)
    label["2_1"] = make_label(name2, name1)
    label["1_1"] = make_label(name1, name1)
    label["2_2"] = make_label(name2, name2)
    for key in score.keys():
        plt.plot(xvals, score[key], clr[key], label=label[key])
    plt.tight_layout()
    plt.ylabel("Score")
    #plt.xlabel("Day Hour in August 2008 (UTC)")
    #dform = DateFormatter("%d %H")
    ax = plt.gca()
    #ax.xaxis.set_major_formatter(dform)
    #ax.set_xticklabels(rotation=45, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    return ax


class FeatureTracker:
    """
    
    """
    def __init__(self, method="bgm", nnn=10, poll=4, lev=[1000, 20000]):
        self.cmap = matplotlib.cm.get_cmap("viridis")
        self.method = method
        self.nnn = nnn
        self.poll = poll
        self.lev = lev

    def plotscore(self):
        runlist = self.runlist
        score = []
        pscore = []
        iii = 0
        for tfit in runlist:
            gfit = tfit.gfit
            score.append(gfit.score(tfit.xra))
            if iii != 0:
                pscore.append(pfit.score(tfit.xra))
            pfit = gfit
            iii += 1
        sns.set_style("whitegrid")
        plt.plot(self.datelist, score, "-k.")
        plt.plot(self.datelist[1:], pscore, "-r.")
        plt.tight_layout()
        return score, pscore

    def create(self, pdump):
        """
        creates a fit for each time period found in the DataFrame.
        """
        runlist, datelist = makelist(pdump, self.method, self.nnn, self.poll, self.lev)
        self.runlist = runlist
        self.datelist = datelist

    def plot1(self, dim="ht", xlim=None, ylim=None):
        Trunlist = self.runlist
        cmap = self.cmap
        iii = 0
        jjj = 0
        mmm = 0
        clr = (0, 1, 0)
        clist = []
        c0done = False
        c0 = 1
        c1 = 0
        nnn = 1
        dd = 1 / 25.0
        nnn = 1 - dd
        for tfit in Trunlist:
            clr2 = cmap(nnn)
            # print(datelist[iii])
            # tfit.scatter(dim=dim,cmap='tab20')
            cnts = tfit.plot_centers(dim=dim, sym="ko", clr=clr2, MarkerSize=4)
            # plt.show()
            ax = plt.gca()
            if dim == "ht":
                ax.set_ylabel("latitude")
                ax.set_xlabel("longitude")
            elif dim == "lat":
                ax.set_ylabel("ht (km)")
                ax.set_xlabel("longitude")
            if xlim:
                ax.set_xlim(xlim[0], xlim[1])
            if ylim:
                ax.set_ylim(ylim[0], ylim[1])
            nnn = nnn - dd
            clr2 = cmap(nnn)
            iii += 1
        plt.savefig("track.png")

    def plot3d(self, xlim=None, ylim=None, zlim=None):
        Trunlist = self.runlist
        datelist = self.datelist
        cmap = self.cmap
        mmm = 0
        datehash = {}
        jjj = 0
        dd = 1 / 25.0
        iii = 1 - dd
        plt.show()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        for tfit in Trunlist:
            # print(datelist[iii])
            # tfit.scatter(dim=dim,cmap='tab20')
            clr2 = cmap(iii)
            datehash[datelist[jjj]] = clr2
            tfit.plot_centers3d(ax, sym="ko", clr=clr2, MarkerSize=20)
            ax = plt.gca()
            ax.set_zlabel("ht (km)")
            ax.set_ylabel("latitude")
            ax.set_xlabel("longitude")
            # plt.show()
            iii = iii - dd
            jjj += 1
            # clr2 = cmap(iii)
        plt.savefig("3dtrack.png")

    def make_legend(self):
        st = []
        d1 = datetime.datetime(2008, 8, 10, 12)
        dt = datetime.timedelta(hours=3)
        while d1 < datetime.datetime(2008, 8, 11, 12):
            keylist.append(d1)
            d1 = d1 + dt

        keylist.append(datetime.datetime(2008, 8, 11, 11))
        # print(keylist)
        make_legend(keylist, datehash)

    def gen_plot(self, dim="lat"):
        Trunlist = self.runlist
        datelist = self.datelist

        dim = "lat"
        if dim == "lat":
            ylim = [0, 18]
            xlim = [-160, -105]
        elif dim == "ht":
            ylim = [30, 80]
            xlim = [-160, -105]
            # plt.figure(figsize=(10,10))
        clrs = ["r", "m", "g", "b", "c", "y", "k"]
        # clrs = ['r','k','g','b']
        sym = ["o", "*", "+", "x", "^"]
        iii = 0
        for tfit in Trunlist:
            dstr = datelist[iii].strftime("%b %d %H:%M UTC")
            self.subplot(tfit, dim, xlim, ylim, dstr)
            iii += 1

    def plot_one(self, iii, dim="ht", xlim=None, ylim=None):
        dstr = self.datelist[iii].strftime("%b %d %H:%M UTC")
        tfit = self.runlist[iii]
        if dim == "lat" and not ylim:
            ylim = [0, 18]
        if dim == "lat" and not xlim:
            xlim = [-160, -105]
        if dim == "ht" and not ylim:
            ylim = [30, 80]
        if dim == "ht" and not xlim:
            xlim = [-160, -105]
        self.subplot(tfit, dim, xlim, ylim, dstr)
        plt.title(dstr)

    def subplot(self, tfit, dim, xlim, ylim, dstr=""):
        if dim == "ht":
            plt.figure(figsize=(6, 6))
        if dim == "lat":
            plt.figure(figsize=(6, 3))
        tfit.scatter(dim=dim, cmap="tab20")
        tfit.plot_centers(dim=dim, MarkerSize=10)
        tfit.plot_gaussians(dim=dim, saturation=0.2)
        ax = plt.gca()
        if dim == "lat":
            ax.set_aspect(1.5, adjustable="box")
            ax.set_ylabel("ht (km)")
            # ax.set_xlabel('longitude')
            ax.text(-155, 16, dstr, bbox=dict(facecolor="white"))

        if dim == "ht":
            ax.set_aspect(1.0, adjustable="box")
            ax.set_ylabel("latitude")
            ax.set_xlabel("longitude")
            ax.text(-155, 75, dstr, bbox=dict(facecolor="white"))

        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlim(xlim[0], xlim[1])
        # plt.tight_layout()
        plt.savefig("tfit{}.{}.png".format(dim, dstr))
