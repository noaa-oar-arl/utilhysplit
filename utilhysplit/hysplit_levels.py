import numpy as np
import matplotlib.pyplot as plt

def findzx(zagl,a,b,c,nlvl):
    dist = b*b-4*a*(c-zagl)
    zx = (-b+dist**0.5)/(2*a)
    zx = np.min([np.max([1,zx]),nlvl])
    return zx


class zindex_example():

    def __init__(self):
        self.zmdl=25000
        self.zz = 0.5
        self.zter = 2000  #Colorado elevation is about 2000 m
        self.nlvl = 71
        self.a = 5
        self.b = 5
        self.c = 0

        self.zagl = (self.zmdl-self.zter)*(1-self.zz)
        self.zaglp = (self.zmdl)*(1-self.zz)

        self.zsg = get_zsg(self.a,self.b,self.c,self.zmdl,self.nlvl)

        zindex = findzx(self.zaglp,self.a,self.b,self.c,self.nlvl)
        print(zindex)

        # need to subtract 1 because index for python arrays starts from 0.
        # while index for fortran arrays starts at 1.
        #hgts = [self.zmdl*(1-x) for x in self.zsg]
        x1 = self.zsg[int(np.floor(zindex))-1]
        x2 = self.zsg[int(np.ceil(zindex))-1]
        target_index = x1 + (x2-x1)*self.zz
        print('target', self.zz, 'value', target_index)
        


    def h_agl(self):
        target = (self.zmdl-self.zter)*(1-self.zz)
        return target

     


def get_zsg(a,b,c,zmdl,nlvl):
    zsg = []
    zter=0
    for nnn in np.arange(1,nlvl+1):
        # zsg values defined for zter=0
        ztemp = 1-(a*nnn*nnn + b*nnn + c)/(zmdl-zter)
        zsg.append(ztemp)
    zsg = np.array(zsg)
    return zsg


def zsg_plots(zmdl, nlvl):
    alist = [30,20,15,10,5,2.5]
    blist = [-25,-15,-5,-10,5,2.5]
    clist = [5,5,0,10,0,0]

    fig = plt.figure(1,figsize=[10,15])
    ax = fig.add_subplot(3,1,1)
    ax2 = fig.add_subplot(3,1,2)
    ax3 = fig.add_subplot(3,1,3)
    
    #sns.set_style('whitegrid')
    sigmalist = []
    xval = range(1,nlvl+1,1)
    for abc in zip(alist,blist,clist):
        sigma = get_zsg(abc[0],abc[1],abc[2],zmdl,nlvl)
        sigmalist.append(sigma)
        lbl = 'a={}, b={}, c={}'.format(abc[0],abc[1],abc[2])
        ax.plot(xval,sigma, marker='.')
        ax2.plot(xval,sigma,marker='.',label=lbl)
        ax3.plot(xval,sigma,marker='.',label=lbl)

    ax.set_ylim(0,1)
    ax2.set_ylim(0.9,1)
    ax2.set_xlim(0,10)
    ax3.set_ylim(0.4,0.6)
    ax3.set_xlim(45,55)
    ax3.set_ylabel('sigma value')
    ax3.set_xlabel('level index')
    #ax3.grid()


    ax.set_xlabel('Level index')
    ax.set_ylabel('sigma value')
    ax2.set_xlabel('Level index')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles,labels)
    plt.title('ZMDL=25000 m')
    plt.show()




def metlvl(plev):

  sfcp = 1013

  # meters per hPa by 100 hPa intervals (100 to 900)
  zph1 = [17.98,14.73,13.09,11.98,11.15,10.52,10.04,9.75,9.88 ]

  # meters per hPa by 10 hPa intervals (10 to 90)
  zph2 = [31.37,27.02,24.59,22.92,21.65,20.66,19.83,19.13,18.51]

  if plev>100:
     iii = int(plev/100.0)-1
     if iii>8: iii=8
     delz = zph1[iii]
  else:
     iii = int(plev/10.0)-1
     if iii>8: iii=8
     delz = zph2[iii]

  zlev = (sfcp-plev)*delz
  return zlev

def est():


  # meters per hPa by 100 hPa intervals (100 to 900)
  zph1 = [17.98,14.73,13.09,11.98,11.15,10.52,10.04,9.75,9.88 ]

  # meters per hPa by 10 hPa intervals (10 to 90)
  zph2 = [31.37,27.02,24.59,22.92,21.65,20.66,19.83,19.13,18.51]

  return zph1, zph2


def hypsometric(p2,p1,tbar):
    rdry = 287.04 #J/Kg-K
    grav = 9.8 #m/s2
    delz = (p2-p1)*rdry*tbar/grav
    return delz

def nasa(p):
    p = p/10 

    #if p>100:
    #    t = (1/5.256)*np.log10(p/101.29)
    #else p<100:
    h = (1.73 - np.log(p/22.65))/0.000157

    return h


def nasalvl(h):
    # h height in meters
    if h<11000:
        T = 15.04 - 0.00649*h
        p = 101.29*((T+273.1)/288.08)**5.256
    elif h<25000:
        T = -56.46
        p = 22.65*np.exp(1.73-0.000157*h)
    else:
        T = -131.21 + 0.00299*h
        p = 2.488*((T+273.1)/216.6)**-11.388
        
    pres = p*10 #convert to hPa
    return pres
