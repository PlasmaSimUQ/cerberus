# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 11:54:08 2011

@author: uqdbond1
"""

import numpy as np

class LandauWave(object):
    
    def __init__(self,rho0, T0, u0, gamma, nx=0, x=[]):
        
        self.rho0 = rho0
        self.T0 = T0
        self.u0 = u0
        self.gamma = gamma
        self.a0 = np.sqrt(gamma)
        
        if not isinstance(x,np.ndarray):
            x = np.array(x)
        
        if nx:
            dx = 2.0/float(nx)
            xmin = -1 + dx/2.0
            xmax = 1.0 - dx/2.0
            self.x = np.linspace(xmin,xmax,nx)
        if x.any():
            self.x = x
            nx = len(x)
        
        self.nx = nx
        
        self.rho = rho0*np.ones((nx))
        self.p = rho0*T0*np.ones((nx))
        self.T = T0*np.ones((nx))
        self.u = u0*np.ones((nx))
        
        self.initialConditions()
        
        return
        
    def initialConditions(self):
        self.u *= np.sin(np.pi*self.x)
        self.rho *= (1.0 + ((self.gamma - 1)*self.u)/(2.0*self.a0))**(2.0/(self.gamma-1))
        self.p *= (1.0 + ((self.gamma - 1)*self.u)/(2.0*self.a0))**((2.0*self.gamma)/(self.gamma-1))
        self.T = self.p/self.rho
    
    def getSol(self, t, steps = 100, epsilon = 1e-12):
        
        #tbreak = np.abs(2/(self.u0*np.pi*(self.gamma + 1)))
        
        #tstop = t_on_tbreak*tbreak
    
        tstop = t
    
        dt = tstop/steps
    
        t = 0.0
        
        for i in range(steps):
            
            t += dt
            
            err = self.getError(t)
            
            count = 0
            while np.all(np.abs(err) > epsilon): 
            
                err = self.getError(t)
                
                errPrime = self.getErrorPrime(t)
                
                self.u -= err/errPrime
                
                count += 1
            
                if count == 10000:
                    print("Landau wave - exit loop, iter >= 10000")
                    break
        
        self.rho = self.rho0*(1.0 + ((self.gamma - 1)*self.u)/(2.0*self.a0))**(2.0/(self.gamma-1))
        self.p = self.rho0*self.T0*(1.0 + ((self.gamma - 1)*self.u)/(2.0*self.a0))**((2.0*self.gamma)/(self.gamma-1))
        self.T = self.p/self.rho
        
        self.t = t
        #self.tbreak = tbreak
        
        
    def getError(self,t):
        arg = self.x - (self.a0 + 0.5*(self.gamma + 1.0)*self.u)*t
        err = self.u0*np.sin(np.pi*arg) - self.u
        
        return err
    
    def getErrorPrime(self,t):
        arg = self.x - (self.a0 + 0.5*(self.gamma + 1.0)*self.u)*t
        errPrime = -1.0 - np.pi*0.5*self.u0*(1.0 + self.gamma)*t*np.cos(np.pi*arg)
        
        return errPrime
    
if __name__ == "__main__":

    import matplotlib.pylab as plt
    
    np.set_printoptions(precision=5, linewidth=250)
    
    x1 = np.linspace(-1.0,1.0,101)
    
    LW1 = LandauWave(1.0,1.0,0.1,5.0/3.0,x=x1)
    tbreak = np.abs(2/(LW1.u0*np.pi*(LW1.gamma + 1)))
    LW1.getSol(0.9*tbreak,steps = 100)
    
    x2 = np.linspace(-.5,0.5,101)
    
    LW2 = LandauWave(1.0,1.0,0.1,5.0/3.0,x=x2)
    LW2.getSol(0.9*tbreak,steps = 100)
    
    plt.plot(LW1.x, LW1.rho)
    plt.plot(LW2.x, LW2.rho,'r')
    
    plt.show()