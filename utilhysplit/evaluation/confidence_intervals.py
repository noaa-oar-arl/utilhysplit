import numpy as np
from scipy.stats import bernoulli
import matplotlib.pyplot as plt
from scipy.stats import norm

from statsmodels.stats.proportion import proportion_confint

"""
2023 Feb 22 
classes and functions related to Bernoulli distributions.
exceedance probabilities computed from ensemble / monte carlo simulation have Bernoulli distribution.


just use the statsmodels packaged
"""


def example(nsize,p):
    # see https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    # methods types are normal, agresti_coull, beta, wilson, jeffreys, binom_test
    methodlist = ['normal','agresti_coull','beta','wilson','jeffreys','binom_test']
    successes = int(p * nsize)
    clist = []
    for method in methodlist:
         ci = proportion_confint(successes,nsize,method=method)
         cilist.append(ci)
    return cilist


def ciplot(nsize, msize):
    nsize = 100


class Bdist:
    def __init__(self,p,nsize,msize=100, alpha=0.05):
        self.nsize = nsize
        self.p = p
        self.msize = msize
        self.alpha = alpha
        self.zalpha = norm.ppf(1-alpha/2.0)

    def __str__(self):
        rstr = ''
        rstr += 'Ensemble sizes {}\n'.format(self.nsize)
        rstr += 'p value {}\n'.format(self.p)
        rstr += 'estimated variance in p {:0.2e} \n'.format(self.rmean.var())
        rstr += 'real variance in p {:0.2e} Standard Deviation  {:0.2e}\n'.format(self.realmeanvar, np.sqrt(self.realmeanvar))
        rstr += 'smallest p value {}\n'.format(1/self.nsize)
        rstr += 'number ensemble members {}\n'.format(self.nsize * self.p)
        return rstr     
   

    def plotA(self):
        fig = plt.figure(figsize=[10,5])
        xxx = np.arange(0,len(self.cilist))
        yyy = self.rmean
        ci = self.cilist
        rci = self.rcilist
        plt.errorbar(xxx,yyy,yerr=rci,fmt='k.',elinewidth=1,capsize=2,ecolor='r')
        plt.errorbar(xxx,yyy,yerr=ci,fmt='k.',elinewidth=1,capsize=2)
        plt.plot([xxx[0],xxx[-1]],[self.p, self.p],'-r') 

    def check_rci(self):
        yes = 0
        no  = 0
        for val in zip(self.cilist, self.rcilist, self.rmean):
            ci = val[0]
            rci = val[1]
            pn = val[2]
            #print(np.abs(self.p-pn),  ci)
            if np.abs(self.p - pn) < rci: 
              yes += 1 
            else:
              no += 1
        print('alpha set at {}'.format(self.alpha))
        print('value within confidence interval {}%  of time'.format(100*yes/(no+yes)))
        return yes, no 


    def check_ci(self, tp='c'):
        yes = 0
        no  = 0
        for val in zip(self.cilist, self.acilist, self.rcilist, self.rmean):
            if tp == 'a':
               ci = val[1]
            elif tp == 'r':
               ci = val[2]
            else:
               ci = val[0]
            pn = val[3]
            #print(np.abs(self.p-pn),  ci)
            if np.abs(self.p - pn) < ci: 
              yes += 1 
            else:
              no += 1
        print('alpha set at {}'.format(self.alpha))
        print('value within confidence interval {}%  of time'.format(100*yes/(no+yes)))
        return yes, no 
 
    def makelists(self): 
        rlist = []
        rmean = []
        rvar = []
        logp = []
        rmeanvar = []
        cilist = []
        acilist = [] # alternative ci list
        rcilist = []
        logcilist = []
        # wilson score interval.
        wsi = []
        wplist = []
        p = self.p
        nsize = self.nsize
        zalpha = self.zalpha        
        for mmm in np.arange(0,self.msize):
            r1 = bernoulli.rvs(p,size=nsize)
            rlist.append(r1)

            # this is the estimate of p
            pn = r1.mean()
            rmean.append(pn)

            # this is the variance of the distribution
            rvar.append(r1.var())

            # confidence interval. 
            # see https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

            ci = zalpha * np.sqrt(pn*(1-pn)/nsize)
            cilist.append(ci)
            if pn==0:
               apn = 1/(nsize+1) 
               ci = zalpha * np.sqrt(apn*(1-apn)/nsize)
            acilist.append(ci)
           
            qqq = zalpha**2/nsize 
            bbb = np.sqrt(pn*(1-pn)/nsize + 0.25*qqq/nsize)
            ci = zalpha/(1+qqq)*bbb
            wsi.append(ci) 
            wp = (pn+0.5*qqq)/(1+qqq)
            wp.append(wp) 
 
            ci = zalpha * np.sqrt(p*(1-p)/nsize)
            rcilist.append(ci)
            

            # this is an estimate of the variance in p.
            rmeanvar.append(pn*(1-pn)/nsize)
 
            #
            logp.append(np.log10(pn))
            ci = zalpha * np.sqrt((1-pn)/(pn*nsize))
            logcilist.append(ci)

        self.rmean = np.array(rmean)
        self.rvar  = np.array(rvar)
        self.realmeanvar = p*(1-p)/nsize
        self.cilist = cilist
        self.rcilist = rcilist
        self.acilist = acilist
        self.logp = logp
        self.logcilist = logcilist

def calcB(p, nsize, msize=100, alpha=0.05):
    nsize = 100
    msize = 100
    rlist = []
    rmean = []
    rvar = []
    rmeanvar = []
    alpha = 0.05
    zalpha = (1-alpha/2.0)
    cilist = []
    for mmm in np.arange(0,msize):
        r1 = bernoulli.rvs(p,size=nsize)
        rlist.append(r1)

        # this is the estimate of p
        pn = r1.mean()
        rmean.append(pn)

        # this is the variance of the distribution
        rvar.append(r1.var())

        # this is the confidence interval from equation 3.
        ci = zalpha * np.sqrt(pn*(1-pn)/nsize)
        cilist.append(ci)

        # this is an estimate of the variance in p.
        rmeanvar.append(pn*(1-pn)/nsize)


    rmean = np.array(rmean)
    rvar  = np.array(rvar)

    # this is the actual variance in p given size nsize.
    realmeanvar = p*(1-p)/nsize
    print('For ensemble sizes {}'.format(nsize))
    print('For p value {}'.format(p))
    print('estimated variance in p {:0.2e}'.format(rmean.var()))
    print('real variance in p {:0.2e}'.format(realmeanvar))
    return 
