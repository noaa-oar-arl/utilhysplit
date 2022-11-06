#!/n-home/alicec/anaconda/bin/python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utilhysplit.plotutils import colormaker
#toy model of dispersion equation used in HYSPLIT.
#given a time step, dt, Lagrangian time scale,tl, and number of iterations.
#starts 10,000 particles from a single place and plots the FWHM of their spatial distribution as a funciton of time.



#def velocity_field(position):
#    return (1,0)


# Notes on HYSPLIT pardsp.f

# HYSPLIT determine vertical time step
# DELZ - distance between vertical levels. 
# wmix is vertical turbulence.
# take max value in the profile
# vscale is the vertical time scale.

# 0.125*DELZ**2 / 
#ZMDL is top of domain.
#ZSG internal model sigma levels.


class SigmaProfile():

    def __init__(self,levs=[0,2000.0,50.0], blh=1500):
        #use a parabola in the boundary layer.
        #AAA = 0.5
        #BBB = 0.5
        AAA = 0.5
        BBB = 0.5
        CCC = 0.05
        space = (levs[1] - levs[0]) / float(levs[2])
        self.zzz = np.arange(levs[0], levs[1], space+1)
        self.sig = []
        self.sigdz = []
        self.compute_sig(AAA,BBB,CCC, blh)

    def add_sig(self, sig,zzz):
        self.sig = sig
        self.zzz = zzz
        self.compute_sigdz()
        return -1 

    def find_sig(self, ht, verbose=False):
        if isinstance(ht, np.ndarray):
           if len(ht) == 1:
               ht = ht[0]
           else:
               print('find_sig input height wrong type', ht)
        #find sig and sigdz based on interpolation between levels.
        sig =self.sig
        sigdz=self.sigdz
        zzz=self.zzz
        zii = np.where(zzz>ht)
        iii = zii[0][0]-1
            
        #print(zii, ht, iii)
        #space = self.zzz[1] - self.zzz[0]
        #iii = int(np.floor(ht / float(space)))
        #iii = np.min([iii, len(zzz)-2])
        if iii < 0: iii=0 
        #if verbose: print('space', space)
        if verbose: print('ht', ht, iii, zzz[iii], zzz[iii+1])
        fact = (zzz[iii+1] - ht) / float(zzz[iii+1]-zzz[iii])
        newsig = fact*(sig[iii+1]-sig[iii]) + sig[iii]
        if verbose: print('zzz', zzz[iii], zzz[iii+1], fact)
        if verbose: print('newsig,sig,sig', newsig, sig[iii], sig[iii+1])
        if verbose: print('sigdz', sigdz[iii], sigdz[iii+1])
        return newsig , sigdz[iii]

    def compute_sig(self, AAA, BBB, CCC, blh):
        """
        creates an artifical sigma prfile using a parabola.
        """

        zzz = self.zzz.copy()
        zzz = zzz /  blh
        sig = -1 * AAA*zzz**2 * zzz + BBB * zzz + CCC
        sig = np.array(sig)*50
 
        vpi = np.where(sig> 0.5)
        sig[vpi] = sig[vpi] * 0.2
         
        vpi = np.where(sig< 0.01)
        sig[vpi] = 0.01
        print('max sig', np.max(sig))
        self.sig = sig
        self.compute_sigdz()
        #vpi = np.where(sig> 0.5)
        #sig[vpi] = sig[vpi] * 0.2
        #plt.plot(sig, zzz2, '--k.')

    def compute_sigdz(self):
        iii=0
        sig = self.sig
        zzz2 = self.zzz
        sigdz = np.array(zzz2) * 0
        for val in sig:
            if iii==0:
               sigdz[iii] = (sig[1]-sig[0])/(zzz2[1]-zzz2[0])
            else:
               sigdz[iii] = (sig[iii]-sig[iii-1])/(zzz2[iii]-zzz2[iii-1])
            iii+=1
        self.sigdz = sigdz
         
    def plot_sig(self, xdz = 100):
        plt.plot(self.sigdz*xdz, self.zzz, '--b.')
        plt.plot(self.sig, self.zzz, '--k.')
        plt.xlabel('Sig (black), sigdz * {} (blue)'.format(xdz))
        #plt.show()
 

def find_sigma(pxra, pyra):
    #make ustar same everywhere to begin with.
    #TO DO vary this with position.
    ustar = 0.01  #m/s  approximate value for ustar
    ##this is Kantha clayson for stble or neutral surface layer
    sigmau = (4.0 * ustar**2) **0.5 * np.ones_like(pxra)
    sigmav = (5.0 * ustar**2) **0.5 * np.ones_like(pxra)
    return  sigmau, sigmav


class Particle():
   
    def __init__(self, xpos, ypos, vel=0, sigw=0.1, sigdz=0):
        self.xpos = xpos           # xpos of particle
        self.ypos = ypos           # ypos of particle
        self.vel = vel             # vel of particle
        self.sigw = sigw
        self.sigdz = sigdz
        self.history = [ypos]      #list of ypositions
        self.time = 0                
        self.time_history = [0]    # list of times
        self.vel_history = [vel]   # list of velocities
        self.sigw_history = [sigw] # list of sigw
        self.sigdz_history = [sigdz] # list of sigdz
        self.delty = [0]             # list of dz (changes in ypos

    def move_simple(self, vely, dt, rcorr):

        BBB=(1-rcorr**2)
        vel2 = rcorr *self.vel + vely *  BBB**0.5
        self.ypos = self.ypos + vel2 * dt
        self.time += dt


        self.time_history.append(self.time)
        self.vel_history.append(vel2)         
        
        self.vel = vel2
        ypos2 = self.ypos + vel2 * dt
        ypos2 = self.ra2float(ypos2)
        dyy = ypos2 - self.ypos
        self.delty.append(dyy)
        if  ypos2 < 0: ypos2= -1 * ypos2

        self.ypos= ypos2
        self.history.append(self.ypos)

    def ra2float(self,val):
        if isinstance(val, np.ndarray):
           if len(val) == 1:
               val = val[0]
           else:
               print('value wrong type ', val)
               val = -1
        return val
 
    def move(self, ww, dt, rcorr, sigw, sigdz, vscale):

        self.vel = self.ra2float(self.vel)
        self.sigw = self.ra2float(self.sigw)
        self.ypos = self.ra2float(self.ypos)

        BBB=(1-rcorr**2)**0.5
        AAA=(1-rcorr)
        wratio = self.vel/self.sigw
        sigw2 = self.sigw + self.vel*dt *self.sigdz
        wratio2 = rcorr*wratio + BBB * ww/sigw2 +  AAA*vscale*sigdz
        #wratio2 = rcorr*wratio + BBB * ww/sigw2
        vel2 = wratio2 * sigw2
        
        self.vel = vel2
        self.vel_history.append(vel2)         
        ypos2 = self.ypos + vel2 * dt
        ypos2 = self.ra2float(ypos2)
        dyy = ypos2 - self.ypos
        self.delty.append(dyy)
        if  ypos2 < 0: ypos2= -1 * ypos2
        #self.delty.append(ypos2-self.ypos)
        self.ypos= ypos2
        self.history.append(self.ypos)
      

        self.time += dt
        self.time_history.append(self.time)

        self.sigdz_history.append(sigdz)
     
        self.sigw = sigw2
        self.sigw_history.append(sigw)

    def plotsigw(self, clr):
        if len(self.sigw_history) != len(self.history):
           print('Cannot plot sigdw vs ypos ')
           return False
        if isinstance(clr,str):
            plt.plot(self.sigw_history, self.history, clr)
        else:
            plt.plot(self.sigw_history, self.history, c = clr, marker='.', linewidth=1)
        plt.xlabel('sigw')   
        plt.ylabel('y position')   
        return True

 
    def plottest(self, clr):
        if len(self.sigdz_history) != len(self.history):
           print('Cannot plot sigdz vs ypos ')
           return False
        if isinstance(clr,str):
            plt.plot(self.sigdz_history, self.history, clr)
        else:
            plt.plot(self.sigdz_history, self.history, c = clr, marker='.', linewidth=1)
        plt.xlabel('sigdz')   
        plt.ylabel('y position')   
        return True

    def plotall(self, clr):
        #print('plotting')
        if len(self.time_history) != len(self.history):
           print('Cannot plot time vs ypos ')
           return False
        if isinstance(clr,str):
            plt.plot(self.time_history, self.history, clr)
        else:
            plt.plot(self.time_history, self.history, c=clr,marker='.', linewidth=0,markersize=1)
        return True

class DispersionSim:

    def __init__(self,numpar=500,dt=1):
        self.numpar=numpar #number of particles to use.
        self.TLu = 5 # seconds vertical lagrangian time scale
        #TLw = 200   # seconds vertical lagrangian time scale
        self.dt = 1    # time step (seconds)
        self.rcorr = np.exp( -1 * float(self.dt) / self.TLu)  #this is the correlation coefficient.
        self.B = (1-self.rcorr**2)                       #this is multiplied by sigma * lambda
        self.sprofile =  SigmaProfile()
        self.sigmax = 0                   #no dispersion in horizontal.
        self.uuu = 0    #mean velocity in u direction everywhere is 2 m/s.
        self.vvv = 0    #mean velocity is v direction everywhere is 0m/s.

        self.parlist = [] # created by create_parlist

    def add_sprofile(self, sig, zzz):
        self.sprofile.add_sig(sig,zzz)



    def reset_dt(self,dt):
        self.dt = dt    # time step (seconds)
        self.rcorr = np.exp( -1 * float(self.dt) / self.TLu)  #this is the correlation coefficient.
        self.B = (1-self.rcorr**2)                       #this is multiplied by sigma * lambda
         

    def create_parlist_uniform(self, ytop=1400):     
        dz = ytop / self.numpar
        parlist = []
        ypos = dz
        for nnn in np.arange(1,self.numpar):
            xpos = nnn
            ypos = dz + nnn*dz
            parlist.append(Particle(xpos, ypos))
        self.parlist = parlist

    
    def create_parlist(self, yinit=1400):     
        ypos = yinit
        parlist = []
        for nnn in np.arange(1,self.numpar):
            xpos = nnn
            parlist.append(Particle(xpos, ypos))
        self.parlist = parlist

    def run(self,iterations,verbose=False):
        #pdist = []
        dplot = []
        for iii in np.arange(0,iterations):
            for part in self.parlist:
                #get the sigw dgdz according to yposition.
                sigw, sigdz = self.sprofile.find_sig(part.ypos,verbose=verbose)
                #pull the zvel from distribution of width sigw
                #print(sigw)
                zvel = np.random.normal(0,np.abs(sigw),1) 
                #move the particle.
                part.move(zvel, self.dt, self.rcorr, sigw, sigdz, self.TLu) 
                #part.move_simple(zvel, self.dt, self.rcorr)
                dplot.extend(part.delty)
        self.dplot = dplot 

    def plot_dplot(self):
        sns.set_style('whitegrid')
        sns.distplot(self.dplot)

    def plotone(self):#fig = plt.Figure(1)
        pdist = []
        print(len(self.parlist))
        if len(self.parlist) <=100:
            cm = colormaker.ColorMaker('viridis', len(self.parlist), ctype='rgb')
            clist = cm() 
        else:
            clist = 'r'
        for iii, part in enumerate(self.parlist):
            if isinstance(clist,list): clr = clist[iii]
            else: clr='r'
            #part.plottest(clr=clr)
            part.plotsigw(clr=clr)
            pdist.append(part.ypos)
        plt.show()

    def plot_all(self):#fig = plt.Figure(1)
        pdist = []
        if len(self.parlist) <=100:
            cm = colormaker.ColorMaker('viridis', len(self.parlist), ctype='rgb')
            clist = cm() 
        else:
            clist = 'r'
        for iii, part in enumerate(self.parlist):
            if isinstance(clist,list): clr = clist[iii]
            else: clr='r'
            part.plotall(clr=clr)
            pdist.append(part.ypos)
        plt.show()
        sns.distplot(pdist,bins=10)
        plt.show()

def dispersion():
    sns.set_style("darkgrid",{"axes.facecolor": "0.8"})
    sns.set_context("talk")
    numpar=500 #number of particles to use.
    TLu = 5 # seconds vertical lagrangian time scale
    #TLw = 200   # seconds vertical lagrangian time scale
    dt = 1    # time step (seconds)
    iterations = 1000 #number of iterations.

    rcorr = np.exp( -1 * float(dt) / TLu)  #this is the correlation coefficient.
    B = (1-rcorr**2)                       #this is multiplied by sigma * lambda
    #print( 'BBBB', B)

    #rand2 = np.random.normal(0,1,nnn) 
    #print rand

    #pxra = np.arange(0,numpar,1)   #spread out particles in the horizontal so can
    #pyra = 100 * np.ones_like(rand) #start all particles at a certain height.
    sprofile =  SigmaProfile()
    sprofile.plot_sig()
    plt.show()
    #sprofile.compute_sig()
    #sprofile.plot_sig()
    sigmax = 0                   #no dispersion in horizontal.

    #Ut = rand2 * sigmau
    #Vt = rand * sigmav
    #Ut=np.zeros_like(rand)  #initial turbulent velocity is zero.
    #Vt=np.zeros_like(rand)  #initial turbulent velocity is zero.
    uuu = 0    #mean velocity in u direction everywhere is 2 m/s.
    vvv = 0    #mean velocity is v direction everywhere is 0m/s.

    #pxra = np.zeros_like(rand)   #initial position of all particles is 0
    #pyra = np.zeros_like(rand)
    #at each time step compute new turbulent velocities and positions of all nnn particles.
    #sigma = find_sigma(pxra, pyra)

    #create list of particle objects. 
    ypos = 1400
    parlist = []
    for nnn in np.arange(1,numpar):
        xpos = nnn
        parlist.append(Particle(xpos, ypos))
    pdist = []
    dplot = []
    for iii in np.arange(0,iterations):
        for part in parlist:
            #get the sigw dgdz according to yposition.
            sigw, sigdz = sprofile.find_sig(part.ypos)
            #pull the zvel from distribution of width sigw
            zvel = np.random.normal(0,sigw,1) 
            #move the particle.
            print('HERE', sigw, sigdz)
            part.move(zvel, dt, rcorr, sigw, sigdz, TLu) 
            dplot.extend(part.delty)
    #fig = plt.Figure(1)
    for part in parlist:
        part.plottest(clr='r.')
        pdist.append(part.ypos)
    plt.show()
    for par in parlist:
        par.plotall(clr='r.')
    plt.show()
    #fig = plt.Figure(2)
    sns.distplot(pdist)
    plt.show()
    sns.distplot(dplot)
    plt.show()
    sys.exit()
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
        print('at iteration {} out of {}'.format(jjj,iterations))
        for iii in np.arange(0,nnn):
            tempu = posu[iii][0:jjj]
            tempv = posv[iii][0:jjj] 
            #print tempu
            #print tempv
            plt.plot(tempu, tempv, c[iii])
        print(jjj)
        ax = plt.gca()
        ax.set_ylim([-1000,1000])
        ax.set_xlim([-10,12000])
        plt.savefig('pos' + str(jjj).zfill(2) + '.jpg')
        plt.show()

#dispersion()
