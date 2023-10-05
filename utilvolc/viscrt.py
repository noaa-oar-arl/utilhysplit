from utilvolc import volcMER
import numpy as np


def get_itab():
    
#                     summit height in kft
#!                 0 02 04 06 08 10 12 14 16 18

    itab = np.array([[15,0,0,0,0,0,0,0,0,0],        # 0
                  [15,15,0,0,0,0,0,0,0,0],          # 2                   
                  [15,15,15,0,0,0,0,0,0,0],         # 4                       
                  [15,15,15,15,0,0,0,0,0,0],        # 6                        
                  [15,15,15,15,15,0,0,0,0,0],       # 8                         
                  [16,15,15,15,15,15,0,0,0,0],      # 10                          
                  [16,16,15,15,15,15,15,0,0,0],     # 12                           
                  [16,16,16,16,15,15,15,15,0,0],    # 14                            
                  [16,16,16,16,16,15,15,15,15,0],   # 16                             
                  [16,16,16,16,16,16,15,15,15,15],  # 18                              
                  [17,16,16,16,16,16,16,15,15,15],  # 20                              
                  [17,17,17,16,16,16,16,16,16,15],  # 22                              
                  [17,17,17,17,17,16,16,16,16,16],  # 24                              
                  [17,17,17,17,17,17,17,16,16,16],  # 26                              
                  [17,17,17,17,17,17,17,17,17,16],  # 28                              
                  [18,18,17,17,17,17,17,17,17,17],  # 30                              
                  [18,18,18,18,18,18,17,17,17,17],  # 32                              
                  [18,18,18,18,18,18,18,18,18,18],  # 34                              
                  [18,18,18,18,18,18,18,18,18,18],  # 36                              
                  [18,18,18,18,18,18,18,18,18,18]]) # 38                            
    return itab


def viscrt(lvlone, lvltwo,ireduc):
    # convert meters to feet.
    ivsh = lvlone*3.2808+0.5 
    iact = lvltwo*3.2808+0.5
    itab = get_itab()
    print(ivsh,iact)

    if iact>=5e4:
       mvis = -19+ireduc
    elif iact>4e4 and iact<5e4:
       mvis=-18+ireduc
    else:
       iii = np.min([ivsh,2e4])/2e3
       jjj = iact/2e4
       iii =int(np.floor(iii))
       jjj =int(np.floor(jjj))
       mvis = -1*itab[iii,jjj] + ireduc
       print('i,j', iii,jjj, mvis)
    return mvis 

def thresh2conc(vent,ht):
    mer = volcMER.mastinMER((ht-vent)/1000.0) #kg/s
    mer2 = mer*3600*1e6 #mg/h
    thresh = viscrt(vent,ht,0)
    print('MER {:1.1e} kg/s  thresh {}'.format(mer,thresh))
    thresh = 10.0**thresh
    return thresh*mer2
  
def check():
    vent=4000
    htlist = np.arange(10000,25000,2000)
    for ht in htlist:
        t = thresh2conc(vent,ht)
        print('vent {}, plume {}, thresh {}'.format(vent/1000.0, ht/1000.0, t))





    

