import numpy as np
import matplotlib.pyplot as plt

class handy(object):

    def __init__(self, deathR_min=1.e-2, deathR_max=7.e-2, 
                       birthR=[3.e-2, 3.e-2], 
                       natRegenR=1.e-2, natCapacity=100.,
                       minSalary=5.e-4, minConsumption=5.e-3, 
                       inequality=1., depletionR=6.67e-6
                       ):
        '''
            s: state variable (X, Y, W)
            X: population (Xw, Xn) => s[0][0]=Xw, s[0][1]=Wn
        '''
    
        self._minConsumption=minConsumption
        self._inequality=inequality
        self._minSalary=minSalary
        self._deathR_min=deathR_min
        self._deathR_max=deathR_max
        self._birthR=birthR
        self._depletionR=depletionR
        self._natRegenR=natRegenR
        self._natCapacity=natCapacity

        self.nClasses=2

        self.isIntegrated=False
    #----| Standard overload |------------------------------

    def __getitem__(self, i):
        if i>=self.nDt+1:
            raise IndexError()
        else:
            return self.x[:,i], self.y[i], self.w[i]

    #----| Equilibrium |------------------------------------



    #def eqDepletionR()


    #----| State functions |--------------------------------

    def minConsumption(self, x, y, w):
        return np.array([self._minConsumption, 
                         self.inequality(x, y, w)*self._minConsumption])
    
    def inequality(self, x, y, w):
        return self._inequality

    def wealthThreshold(self, x, y, w):
        return (self.minConsumption(x, y, w)[0]*x[0]
                    +self.minConsumption(x, y, w)[1]*x[1])

    def salary(self, x, y, w):
        return np.array([self._minSalary, 
                         self.inequality(x, y, w)*self._minSalary])

    def consumption(self, x, y, w):
        return np.array([np.min([1., w/self.wealthThreshold(x, y, w)])
                                *self.salary(x, y, w)[0]*x[0],
                         np.min([1., w/self.wealthThreshold(x, y, w)])
                                *self.salary(x, y, w)[1]*x[1]])


    def deathR(self, x, y, w):
        famineFactor=np.zeros(self.nClasses)
        for i in xrange(self.nClasses):
            if x[i]==0:
                famineFactor[i]=1.
            else:
                famineFactor[i]=(1.
                    -self.consumption(x,y,w)[i]/(self._minSalary*x[i]))

        return self._deathR_min+(np.array([np.max([0., famineFactor[0]]),
                                         np.max([0., famineFactor[1]])])
                                *(self._deathR_max-self._deathR_min))

    def birthR(self, x, y, z):
        return np.array(self._birthR)
        
    def depletionR(self, x, y, w):
        return self._depletionR

    def workProd(self, x, y, w):
        return self._depletionR*x[0]*y

    def natureProd(self, x, y, w):
        return self._natRegenR*y*(self._natCapacity-y)


    #----| Integration |------------------------------------


    def increment(self, x, y, w):
        dx=(self.birthR(x,y,w)-self.deathR(x, y, w))*x
        dy=self.natureProd(x,y,w)-self.workProd(x,y,w)
        dw=self.workProd(x,y,w)-self.consumption(x,y,w).sum()
        return dx, dy, dw

    def initialize(self, x0, y0, w0, tInt, dt=1.):
        self.nDt=int(tInt/dt)
        
        # state variables
        self.x=np.zeros(shape=(self.nClasses, self.nDt+1,))
        self.y=np.zeros(self.nDt+1)
        self.w=np.zeros(self.nDt+1)
        self.time=np.zeros(self.nDt+1)
        self.x[:,0]=x0
        if y0==None: self.y[0]=self._natCapacity
        else: self.y[0]=y0
        self.w[0]=w0

        self.time[0]=0.

    def integrate(self, tInt, dt=1., x0=[100.,0.], y0=None, w0=0.):
        

        # initialization
        self.initialize(x0, y0, w0, tInt, dt=dt)

        for i in xrange(self.nDt):
            dx, dy, dw=self.increment(self.x[:,i], 
                                      self.y[i],
                                      self.w[i])
            self.x[:,i+1]=self.x[:,i]+dx*dt
            self.y[i+1]=self.y[i]+dy*dt
            self.w[i+1]=self.w[i]+dw*dt

            self.time[i+1]=(i+1)*dt
        self.isIntegrated=True

            
    #----| Post-processing |--------------------------------

    def plotPop(self, axe=None, showEQ=False):
        if not self.isIntegrated: raise Exception('Not integrated')

        color_l=['b', 'r', 'm']
        if axe==None:
            axe=plt.subplot(111)
        for i in xrange(self.nClasses):
            axe.plot(self.time, self.x[i], 
                     color=color_l[np.mod(i,len(color_l))],
                     label='$x_{%d}$'%i)
            if i>0 and showEQ:
                axe.plot(self.time, 
                        self.inequality(self.x[i], 
                                        self.y[i], 
                                        self.w[i])*self.x[i], 
                        color=color_l[np.mod(i,len(color_l))],
                        linestyle='--',
                        label='$x_{%d}^{EQ}$'%i)
        axe.set_xlabel('$t$')
        axe.set_ylabel('population')
        axe.legend(loc='upper left')
        return axe
        

    def plotDemo(self, axe=None):
        if not self.isIntegrated: raise Exception('Not integrated')
        axe=self.plotPop(axe=axe, showEQ=False)
        axe2=axe.twinx()
        deathR=np.zeros(shape=(self.nClasses, self.nDt+1))
        birthR=np.zeros(shape=(self.nClasses, self.nDt+1))
       
        for i in xrange(self.nDt+1):
            deathR[:,i]=self.deathR(self.x[:,i], self.y[i], self.w[i])
            birthR[:,i]=self.birthR(self.x[:,i], self.y[i], self.w[i])
        
        color_l=['b', 'r', 'm']
        for i in xrange(self.nClasses):
           axe2.plot(self.time, birthR[i], 
                     color=color_l[np.mod(i,len(color_l))],
                     linestyle=':',
                     label=r'$\beta_{%d}$'%i)
           axe2.plot(self.time, deathR[i], 
                     color=color_l[np.mod(i,len(color_l))],
                     linestyle='--',
                     label=r'$\alpha_{%d}$'%i)
        axe2.axhline(y=self._deathR_min, label=r'$\alpha_m$',
                     color='k', linestyle='-')
        axe2.axhline(y=self._deathR_max, label=r'$\alpha_M$',
                     color='k', linestyle='--')
        axe2.set_ylabel('$y^{-1}$')
        axe2.legend(loc='upper right')

        return axe, axe2
        

    def plotAll(self, axe=None, showEQ=False):
        if not self.isIntegrated: raise Exception('Not integrated')
        axe=self.plotPop(axe=axe, showEQ=showEQ)
        axe2=axe.twinx()
        axe2.plot(self.time, self.y, label='$y$', color='g')
        axe2.plot(self.time, self.w, label='$w$', color='c')
        axe2.set_ylabel('eco$')
        axe2.legend(loc='upper right')

        return axe, axe2
