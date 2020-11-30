#!/n-home/alicec/anaconda/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
import seaborn as sns
import sys
from sklearn.mixture import GaussianMixture as GMM
from sklearn.cluster import KMeans
from utilhysplit import par2conc
#toy model of dispersion equation used in HYSPLIT.
#given a time step, dt, Lagrangian time scale,tl, and number of iterations.
#starts 10,000 particles from a single place and plots the FWHM of their spatial distribution as a funciton of time.



#def velocity_field(position):
#    return (1,0)

#def apply_kmeans(nnn,xra):
#    kpredict = KMeans(n_clusters=nnn,random_state=0).fit_predict(xra)
#    plt.scatter

def add_bandwidth(xra,clr):
    patches = []
    for pos in xra:
        patches.append(Ellipse(pos,width=400,height=100,fill=True,color=clr,alpha=0.3))
    pc = PatchCollection(patches,alpha=0.3)
    return pc

def apply_gmm(nnn,xra,figname='gmm.png'):
    # xra should be numpy array list of (xpos,ypos) points
    gmm = par2conc.get_gmm(nnn)
    mfit = par2conc.MassFit(gmm,xra) 
    return mfit
    #mfit.scatter(dim='ht')
    #mfit.plot_gaussians(dim='ht')
    #plt.savefig(figname)

def make_poster_movie():
    c1 = sns.xkcd_rgb['royal blue']
    #c1 = sns.xkcd_rgb['dark pink']
    c2 = sns.xkcd_rgb['hot pink']
    numiter=50
    td1 = ToyDispersion()
    td1.nnn=100
    td1.iterations=numiter
    iiib=0
    title = '100 particles. Each particle has 10g mass'
    td1.calculate()
    iii = td1.make_plots(ttl=title,clr=c1)
    td2 = ToyDispersion()
    td2.nnn=1000
    td2.iterations=numiter
    #iiib=numiter
    title = '1000 particles. Each particle has 1g mass'
    td2.calculate()
    print('starting at {}'.format(iii))
    td2.make_plots(iiib=iii,ttl=title,clr=c2)




def find_sigma(pxra, pyra):
    #make ustar same everywhere to begin with.
    #TO DO vary this with position.
    ustar = 0.01  #m/s  approximate value for ustar
    ustar = 0.1  #m/s  approximate value for ustar
    ##this is Kantha clayson equation for stable or neutral surface layer
    sigmau = (4.0 * ustar**2) **0.5 * np.ones_like(pxra)
    sigmav = (5.0 * ustar**2) **0.5 * np.ones_like(pxra)
    return  sigmau, sigmav


class ToyDispersion:

    def __init__(self):
        self.nnn=200 #number of particles to use.
        self.TLu = 10800 # seconds Horizontal lagrangian time scale
        #TLw = 200   # seconds vertical lagrangian time scale
        self.dt = 1*60    # time step (x minutes * 60 seconds/minute)
        self.iterations = 50 #number of iterations.
        self.uuu = 2    #mean velocity in u direction everywhere is 2 m/s.
        self.vvv = 0    #mean velocity is v direction everywhere is 0m/s.



    def calculate(self):
        self.get_corr()
        self.get_rand()
        self.fillra()

    def get_rand(self):
        self.rand = np.random.normal(0,1,self.nnn) 
        self.rand2 = np.random.normal(0,1,self.nnn) 
        self.pxra = np.zeros_like(self.rand)   #initial position of all particles is 0
        self.pyra = np.zeros_like(self.rand)

    def get_corr(self):
        self.rcorr = np.exp( -1 * float(self.dt) / self.TLu)  #this is the correlation coefficient.
        self.B = (1-self.rcorr**2)                       #this is multiplied by sigma * lambda

    def init_hash(self):
        vel={}
        posu={}
        posv={}
        for jjj in np.arange(0,self.nnn):
            vel[jjj] = []
            posu[jjj]=[]
            posv[jjj]=[]
        return vel, posu, posv

    def sub_plot(self,ax, xxx,yyy,clr,ttl):
        plt.plot(xxx, yyy, clr, marker='.',ls='') 
        ax.set_ylim([-1000,1000])
        ax.set_xlim([-10,8000])
        ax.set_xlabel('x position (arbitrary unit)')
        ax.set_ylabel('y position (arbitrary unit)')
        plt.title(ttl,size=12)
        #plt.savefig('traj' + str(iii).zfill(2) + '.jpg') 
        #plt.savefig('ldisp' + str(iii).zfill(2) + '.jpg') 
        plt.tight_layout()

    def fit_gmm(self,nnn,jjj,ht):
        xra = np.array(self.make_xra(jjj,ht))
        mfit = apply_gmm(nnn,xra)
        return mfit

    def make_xra(self,jjj,ht):
        pxra = [self.posu[x] for x in self.posu.keys()]
        pyra = [self.posv[x] for x in self.posv.keys()]
        xxx = [x[int(jjj)] for x in pxra] 
        yyy = [x[int(jjj)] for x in pyra]
        htra = np.ones_like(np.array(xxx))
        htra = htra * ht
        return list(zip(xxx,yyy,htra)) 

    def special_plot(self,ax,m1=1000,m2=250):
        plt.grid(True, which='both') 
        txt1 = 'Histogram method \nto find mass distribution\n'
        txt2 = 'Count particles \n in each grid box, \n divide total mass by volume'
        ax.text(0.05,0.7,txt1+txt2,transform=ax.transAxes,size=10,
               bbox=dict(facecolor='white'))
        ax.xaxis.set_minor_locator(MultipleLocator(m1))
        ax.yaxis.set_minor_locator(MultipleLocator(m2))

    def make_gmm_plotB(self,nnn,iii,sp,ht,ttl='',clr='b',fignum=1,figname='toy_gmm'):
        sns.set_style("darkgrid",{"axes.facecolor": "0.8"})
        sns.set_context("paper")
        fig = plt.figure(fignum)
        ax = plt.gca()
        pxra = [self.posu[x] for x in self.posu.keys()]
        pyra = [self.posv[x] for x in self.posv.keys()]
        xxx = [x[int(iii)] for x in pxra] 
        yyy = [x[int(iii)] for x in pyra]
        xxx2 = [x[int(iii-sp)] for x in pxra] 
        yyy2 = [x[int(iii-sp)] for x in pyra]
        xxx.extend(xxx2)
        yyy.extend(yyy2)
        xra = np.array(list(zip(xxx,yyy)))
        mfit = apply_gmm(nnn,xra)
        #mfit = self.fit_gmm(nnn,iii,ht)
        self.sub_plot(ax,xxx,yyy,clr,ttl)
        mfit.plot_gaussians(dim='ht')
        plt.ylim(-2000,2000)
        plt.grid(False, which='both') 
        #ax.text(0.05,0.1,'Gaussian Mixture Model ',
        #       transform=ax.transAxes,size=10,
        #       bbox=dict(facecolor='white'))
        plt.savefig('{}.png'.format(figname))

    def make_gmm_plot(self,nnn,iii,ht,ttl='',clr='b',fignum=1,figname='toy_gmm'):
        sns.set_style("darkgrid",{"axes.facecolor": "0.8"})
        sns.set_context("paper")
        fig = plt.figure(fignum)
        ax = plt.gca()
        pxra = [self.posu[x] for x in self.posu.keys()]
        pyra = [self.posv[x] for x in self.posv.keys()]
        xxx = [x[int(iii)] for x in pxra] 
        yyy = [x[int(iii)] for x in pyra]
        mfit = self.fit_gmm(nnn,iii,ht)
        self.sub_plot(ax,xxx,yyy,clr,ttl)
        mfit.plot_gaussians(dim='ht')
        plt.ylim(-2000,2000)
        plt.grid(False, which='both') 
        #ax.text(0.05,0.1,'Gaussian Mixture Model ',
        #       transform=ax.transAxes,size=10,
        #       bbox=dict(facecolor='white'))
        plt.savefig('{}.png'.format(figname))

    def make_plots(self,tag='pos',iiib=0,ttl='',clr='b',notshow=True):
        pxra = [self.posu[x] for x in self.posu.keys()]
        pyra = [self.posv[x] for x in self.posv.keys()]
        special = iiib + self.iterations - 1
        jjj=0
        fignum=1
        repeat = 20
        for iii in np.arange(iiib,iiib+self.iterations):
            sns.set_style("darkgrid",{"axes.facecolor": "0.8"})
            sns.set_context("paper")
            fig = plt.figure(fignum)
            xxx = [x[int(jjj)] for x in pxra] 
            yyy = [x[int(jjj)] for x in pyra]
            ax = plt.gca()
            self.sub_plot(ax,xxx,yyy,clr,ttl)
            if iii == special:
               self.special_plot(ax)
            else:
               plt.grid(None)
            plt.savefig('{}.{:03d}.png'.format(tag,iii))     
            if notshow: plt.close()
            fignum += 1 
            jjj+=1
        # repeat the last figure.
        for nnn in np.arange(iii,iii+repeat):
            print('n repeat at {}'.format(nnn))
            fig = plt.figure(fignum)
            ax = plt.gca()
            self.sub_plot(ax,xxx,yyy,clr,ttl)
            self.special_plot(ax)
            plt.savefig('{}.{:03d}.png'.format(tag,nnn))     
            plt.close
        for iii in np.arange(nnn,nnn+repeat):
            print('i repeat at {}'.format(iii))
            fig = plt.figure(fignum)
            ax = plt.gca()
            self.sub_plot(ax,xxx,yyy,clr,ttl)
            self.special_plot(ax,m1=500,m2=125)
            ax.text(0.05,0.1,'Higher spatial resolution',
               transform=ax.transAxes,size=10,
               bbox=dict(facecolor='white'))
            plt.savefig('{}.{:03d}.png'.format(tag,iii))     
            plt.close()
        for nnn in np.arange(iii,iii+repeat):
            print('b repeat at {}'.format(iii))
            fig = plt.figure(fignum)
            ax = plt.gca()
            self.sub_plot(ax,xxx,yyy,clr,ttl)
            self.special_plot(ax,m1=500,m2=125)
            ax.text(0.05,0.1,'A KDE adds a mass distribution \nonto each particle',
               transform=ax.transAxes,size=10,
               bbox=dict(facecolor='white'))
            pc = add_bandwidth(zip(xxx,yyy),'b')
            ax.add_collection(pc)
            plt.savefig('{}.{:03d}.png'.format(tag,nnn))     
            plt.close()
        return nnn
 
    def fillra(self,verbose=False):
        # posu x position for each particle
        # posv v position for each particle
        # vel turbulent velocity component in U direction. 
        vel,posu,posv = self.init_hash()
        pxra = np.zeros_like(self.rand)   #initial position of all particles is 0
        pyra = np.zeros_like(self.rand)
        for iii in np.arange(0,self.iterations):
            rand = np.random.normal(0,1,self.nnn) 
            rand2 = np.random.normal(0,1,self.nnn) 
            # Find the widths of the Gaussian Velocity Distribution. 
            sigmau, sigmav = find_sigma(pxra, pyra)
            Ut = self.rand2 * sigmau
            Vt = self.rand * sigmav
            Vt2 = self.rcorr*Vt + sigmav* self.rand * self.B**0.5  #this is the turbulent velocity component in V direction.
            Ut2 = self.rcorr*Ut + sigmau* self.rand2 * self.B**0.5  #this is the turbulent velocity component in U direction.
            jjj=0
            if verbose:
                print(iii)
                print(pxra)
                print('------')
                print(pyra)
                print('**------')
                print(Ut2)
                print('********************')
            for ui in Ut2:
                vel[jjj].append(ui)
                posv[jjj].append(pyra[jjj])
                posu[jjj].append(pxra[jjj])
                jjj+=1
            # very simple way to calculate new position.
            pyra = pyra + (Vt2+self.vvv) * self.dt
            pxra = pxra + (Ut2+self.uuu) * self.dt
            Vt = Vt2
            Ut = Ut2
        self.vel=vel
        self.posu = posu
        self.posv = posv 



def dispselersion():
    sns.set_style("darkgrid",{"axes.facecolor": "0.8"})
    sns.set_context("talk")
    nnn=200 #number of particles to use.
    TLu = 10800 # seconds Horizontal lagrangian time scale
    #TLw = 200   # seconds vertical lagrangian time scale
    dt = 1*60    # time step (x minutes * 60 seconds/minute)
    iterations = 50 #number of iterations.

    rcorr = np.exp( -1 * float(dt) / TLu)  #this is the correlation coefficient.
    B = (1-rcorr**2)                       #this is multiplied by sigma * lambda
    rand = np.random.normal(0,1,nnn) 
    rand2 = np.random.normal(0,1,nnn) 
    #print rand

    pxra = np.zeros_like(rand)   #initial position of all particles is 0
    pyra = np.zeros_like(rand)

    sigmau, sigmav = find_sigma(pxra, pyra)
    Ut = rand2 * sigmau
    Vt = rand * sigmav
    #Ut=np.zeros_like(rand)  #initial turbulent velocity is zero.
    #Vt=np.zeros_like(rand)  #initial turbulent velocity is zero.
    uuu = 2    #mean velocity in u direction everywhere is 2 m/s.
    vvv = 0    #mean velocity is v direction everywhere is 0m/s.

    #pxra = np.zeros_like(rand)   #initial position of all particles is 0
    #pyra = np.zeros_like(rand)
    #at each time step compute new turbulent velocities and positions of all nnn particles.
    #sigma = find_sigma(pxra, pyra)
   
    vel={}
    posu={}
    posv={}
    for jjj in np.arange(0,nnn):
        vel[jjj] = []
        posu[jjj]=[]
        posv[jjj]=[]

    for iii in np.arange(0,iterations):
        plt.plot(pxra, pyra, '.') 
        rand = np.random.normal(0,1,nnn) 
        rand2 = np.random.normal(0,1,nnn) 
        # Find the widths of the Gaussian Velocity Distribution. 
        sigmau, sigmav = find_sigma(pxra, pyra)
        Vt2 = rcorr*Vt + sigmav* rand * B**0.5  #this is the turbulent velocity component in V direction.
        Ut2 = rcorr*Ut + sigmau* rand2 * B**0.5  #this is the turbulent velocity component in U direction.
        jjj=0
        for ui in Ut2:
            vel[jjj].append(ui)
            posv[jjj].append(pyra[jjj])
            posu[jjj].append(pxra[jjj])
            jjj+=1
        
        #Ut2=0  
        #Vt2=0  
      #iii += 1
        # very simple way to calculate new position.
        pyra = pyra + (Vt2+vvv) * dt
        pxra = pxra + (Ut2+uuu) * dt
        Vt = Vt2
        Ut = Ut2

        #plt.plot(pxra, pyra, 'k.') 
        ax = plt.gca()
        ax.set_ylim([-1000,1000])
        ax.set_xlim([-10,8000])
        ax.set_xlabel('x position (arbitrary unit)')
        ax.set_ylabel('y position (arbitrary unit)')
        #plt.savefig('traj' + str(iii).zfill(2) + '.jpg') 
        plt.savefig('ldisp' + str(iii).zfill(2) + '.jpg') 
        plt.clf()
        #plt.show() 

    for jjj in np.arange(0,nnn):
        plt.plot(vel[jjj])
        plt.savefig('velocity_traces.jpg')
    plt.show()

    c = np.zeros(10)
    c= c.tolist()
    c[1] = sns.xkcd_rgb["azure"]
    c[2] = sns.xkcd_rgb["royal"]
    c[3] = sns.xkcd_rgb["bright blue"]
    c[4] = sns.xkcd_rgb["sky blue"]
    c[5] = sns.xkcd_rgb["denim blue"]
    c[6] = sns.xkcd_rgb["cobalt"]
    c[7] = sns.xkcd_rgb["cornflower"]
    c[8] = sns.xkcd_rgb["marine blue"]
    c[9] = sns.xkcd_rgb["prussian blue"]
    c[0] = sns.xkcd_rgb["purpleish blue"]

    for jjj in  np.arange(0,iterations):
        for iii in np.arange(0,nnn):
            tempu = posu[iii][0:jjj]
            tempv = posv[iii][0:jjj] 
            #print tempu
            #print tempv
            plt.plot(tempu, tempv, c[iii])
        #print jjj
        ax = plt.gca()
        ax.set_ylim([-1000,1000])
        ax.set_xlim([-10,12000])
        plt.savefig('pos' + str(jjj).zfill(2) + '.jpg')
        #plt.show()

#dispersion = ToyDispersion()
#dispersion.calculate()
#dispersion.make_plots()

#dispersion()
#make_poster_movie()
