import numpy as np
import matplotlib.pyplot as plt
from handy import handy


def checkPlot(model, axe):
    if not model.isIntegrated: raise Exception('Not integrated')
    if axe==None:
        axe=plt.subplot(111)
    return axe
    

def plotPop(model, axe=None, showEQ=False):
    axe=checkPlot(model, axe)

    color_l=['b', 'r', 'm']
    for i in xrange(model.nClasses):
        axe.plot(model.time, model.x[i], 
                 color=color_l[np.mod(i,len(color_l))],
                 label='$x_{%d}$'%i)
        if i>0 and showEQ:
            axe.plot(model.time, 
                    model.inequality(model.x[i], 
                                    model.y[i], 
                                    model.w[i])*model.x[i], 
                    color=color_l[np.mod(i,len(color_l))],
                    linestyle='--',
                    label='$x_{%d}^{EQ}$'%i)
    axe.set_xlabel('$t$')
    axe.set_ylabel('population')
    return axe
    

def plotConsumption(model, axe=None):
    axe=checkPlot(model, axe)

    consumption=np.zeros(shape=(model.nClasses, model.nDt+1))    
    for i in xrange(model.nDt+1):
        consumption[:,i]=model.consumption(model.x[:,i], 
                                          model.y[i], 
                                          model.w[i])
    
    color_l=['b', 'r', 'm']
    for i in xrange(model.nClasses):
       axe.plot(model.time, consumption[i], 
                 color=color_l[np.mod(i,len(color_l))],
                 linestyle='-',
                 label=r'$C_{%d}$'%i)
    axe.set_ylabel(r'eco\$/$y$')
    return axe
    

def plotBDRate(model, axe=None):
    axe=checkPlot(model, axe)

    deathR=np.zeros(shape=(model.nClasses, model.nDt+1))
    birthR=np.zeros(shape=(model.nClasses, model.nDt+1))
    for i in xrange(model.nDt+1):
        deathR[:,i]=model.deathR(model.x[:,i], model.y[i], model.w[i])*100.
        birthR[:,i]=model.birthR(model.x[:,i], model.y[i], model.w[i])*100.
    
    color_l=['b', 'r', 'm']
    for i in xrange(model.nClasses):
       axe.plot(model.time, birthR[i], 
                 color=color_l[np.mod(i,len(color_l))],
                 linestyle=':',
                 label=r'$\beta_{%d}$'%i)
       axe.plot(model.time, deathR[i], 
                 color=color_l[np.mod(i,len(color_l))],
                 linestyle='--',
                 label=r'$\alpha_{%d}$'%i)
    axe.axhline(y=model._deathR_min*100., label=r'$\alpha_m$',
                 color='k', linestyle='-')
    axe.axhline(y=model._deathR_max*100., label=r'$\alpha_M$',
                 color='k', linestyle='--')
    axe.set_ylabel('$\%/y$')
    return axe
    
def plotRessources(model, axe=None):
    axe=checkPlot(model, axe)
    axe.plot(model.time, model.y, label='$y$', color='g')
    axe.plot(model.time, model.w, label='$w$', color='c')
    axe.set_ylabel('eco$')
    return axe

def plotState(model, axe=None, showEQ=False):
    axe=checkPlot(model, axe)
    axe=plotPop(model, axe=axe, showEQ=showEQ)
    ax2=axe.twinx()
    ax2=plotRessources(model)
    return axe, ax2
