# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 12:45:02 2011

@author: uqdbond1
"""

# riemann_solver.py

import numpy as np
from math import sqrt
from pylab import *

class riemannSolver(object):
    """
    solves the 1D shock tube Riemann problem
    """
    
    def __init__(self, x, diaphragm_position, gamma, R, 
                 D_L, U_L, P_L, D_R, U_R, P_R, U_scale, D_scale, P_scale, T_scale):
    
    
        # DOMLEN   : Domain length
        # DIAPH1   : Position of diaphragm 1
        # CELLS    : Number of computing cells
        # GAMMA    : Ratio of specific heats
        # TIMEOU   : Output time
        # DL       : Initial density  on left state
        # UL       : Initial velocity on left state
        # PL       : Initial pressure on left state
        # DR       : Initial density  on right state
        # UR       : Initial velocity on right state
        # PR       : Initial pressure on right state
        # USCALE   : Normalising constant
        # DSCALE   : Normalising constant
        # PSCALE   : Normalising constant
        # TSCALE   : Normalising constant
        # R        : Gas Constant
            
        
        #initialise

        

        self.DL = D_L
        self.UL = U_L
        self.PL = P_L
        self.CL = 0.0
        self.DR = D_R
        self.UR = U_R
        self.PR = P_R
        self.CR = 0.0
        self.GAMMA = gamma
        self.R = R
        
        self.X = x
        self.DOMLEN = x[-1] - x[0]
        self.DIAPH1 = diaphragm_position
        self.CELLS = len(x)
        self.GAMMA = gamma
        
        self.USCALE = U_scale
        self.DSCALE = D_scale
        self.PSCALE = P_scale
        self.TSCALE = T_scale
        
        self.G1 = (gamma - 1.0)/(2.0*gamma)
        self.G2 = (gamma + 1.0)/(2.0*gamma)
        self.G3 = 2.0*gamma/(gamma - 1.0)
        self.G4 = 2.0/(gamma - 1.0)
        self.G5 = 2.0/(gamma + 1.0)
        self.G6 = (gamma - 1.0)/(gamma + 1.0)
        self.G7 = (gamma - 1.0)/2.0
        self.G8 = gamma - 1.0
        
       
    def compute(self, timeout):
        """
        calculate the profile for the time specified
        """
        
        self.TIMEOU = timeout
        
        DS = 0.0
        PM = 0.0
        PS = 0.0
        UM = 0.0
        US = 0.0
        
         #Compute sound speeds

        self.CL = sqrt(self.GAMMA*self.PL/self.DL)
        self.CR = sqrt(self.GAMMA*self.PR/self.DR)
        
        #The pressure positivity condition is tested for
        
        if(self.G4*(self.CL+self.CR)<=(self.UR-self.UL)):
            return
            
        
        # Exact solution for pressure and velocity in star
        # region is found
        
        PM, UM = self.STARPU(PM)
        
        #Complete solution at time TIMEOU is found
        
        DATA = np.zeros((self.CELLS,5))

        for I, XPOS in enumerate(self.X):
            S    = (XPOS - self.DIAPH1)/self.TIMEOU
            
            #Solution at point (X,T) = ( XPOS - DIAPH1,TIMEOU) is found
            PM, UM, S, DS, US, PS = self.SAMPLE(PM, UM, S)
            
            #Exact solution profiles are written to data.
            DATA[I,0] = XPOS
            DATA[I,1] = DS/self.DSCALE
            DATA[I,2] = US/self.USCALE
            DATA[I,3] = PS/self.PSCALE
            DATA[I,4] = PS/(DS*self.R)/self.TSCALE
            
        
        return DATA
    
    def STARPU(self,P):
        """
        Purpose: to compute the solution for pressure and
        velocity in the Star Region
        
        global DL UL PL CL DR UR PR CR
        
        """
        
        #initialise values
        FL = 0.0
        FR = 0.0
        FLD = 0.0
        FRD = 0.0
        
        TOLPRE = 1e-6
        
        # Guessed value PSTART is computed
        
        PSTART = self.GUESSP()
        
        POLD  = PSTART
        UDIFF = self.UR - self.UL
        
        for I in range(10):
            FL, FLD = self.PREFUN(POLD, self.DL, self.PL, self.CL)
            FR, FRD = self.PREFUN(POLD, self.DR, self.PR, self.CR)
            P = POLD - (FL + FR + UDIFF)/(FLD + FRD)
            CHANGE = 2.0*abs((P - POLD)/(P + POLD))
            
            if CHANGE <= TOLPRE:
                break
            elif P < 0:
                P = TOLPRE
            POLD = P
        
        # Compute velocity in Star Region
        
        U = 0.5*(self.UL + self.UR + FR - FL)
        
        return P, U

    def GUESSP(self):
        """
        %      Purpose: to provide a guessed value for pressure
        %               PM in the Star Region. The choice is made
        %               according to adaptive Riemann solver using
        %               the PVRS, TRRS and TSRS approximate
        %               Riemann solvers. See Sect. 9.5 of Chapt. 9
        %               of Ref. 1
        
        global G1 G3 G4 G5 G6 G7
        global DL UL PL CL DR UR PR CR

        """
        
        QUSER = 2.0
        
        # Compute guess pressure from PVRS Riemann solver
        
        CUP  = 0.25*(self.DL + self.DR)*(self.CL + self.CR)
        PPV  = 0.5*(self.PL + self.PR) + 0.5*(self.UL - self.UR)*CUP
        PPV  = max([0.0, PPV])
        PMIN = min([self.PL,  self.PR])
        PMAX = max([self.PL,  self.PR])
        QMAX = PMAX/PMIN
        
        if (QMAX <= QUSER) & (PMIN <= PPV) & (PPV <= PMAX):
            PM = PPV
        else:
            if PPV < PMIN:
                #Select Two-Rarefaction Riemann solver
                PQ  = (self.PL/self.PR)**self.G1;
                UM  = (PQ*self.UL/self.CL + self.UR/self.CR + self.G4*(PQ - 1.0))/(PQ/self.CL + 1.0/self.CR)
                PTL = 1.0 + self.G7*(self.UL - UM)/self.CL
                PTR = 1.0 + self.G7*(UM - self.UR)/self.CR
                PM  = 0.5*(self.PL*PTL**self.G3 + self.PR*PTR**self.G3)
            else:
                #Select Two-Shock Riemann solver with PVRS as estimate
                GEL = sqrt((self.G5/self.DL)/(self.G6*self.PL + PPV))
                GER = sqrt((self.G5/self.DR)/(self.G6*self.PR + PPV))
                PM  = (GEL*self.PL + GER*self.PR - (self.UR - self.UL))/(GEL + GER)
        
        return PM
    
    def PREFUN(self,P,DK,PK,CK):
        """
        %      Purpose: to evaluate the pressure functions
        %               FL and FR in exact Riemann solver
        %               and their first derivatives
        
        global G1 G2 G4 G5 G6
        """
        
        if P <= PK:
            #Rarefaction wave
            PRATIO = P/PK;
            F    = self.G4*CK*(PRATIO**self.G1 - 1.0)
            FD   = (1.0/(DK*CK))*PRATIO**(-self.G2)
        else:
            #shock wave
            AK  = self.G5/DK
            BK  = self.G6*PK
            QRT = sqrt(AK/(BK + P))
            F   = (P - PK)*QRT
            FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT
    
        return F,FD
    
    def SAMPLE(self, PM, UM, S):
        """
        % Purpose: to sample the solution throughout the wave
        %               pattern. Pressure PM and velocity UM in the
        %               Star Region are known. Sampling is performed
        %               in terms of the 'speed' S = X/T. Sampled
        %               values are D, U, P
        
        global GAMMA G1 G2 G3 G4 G5 G6 G7
        global DL UL PL CL DR UR PR CR
        """
        
        if S <= UM:
            #Sampling point lies to the left of the contact discontinuity
            if PM <= self.PL:
                #left rarefaction
                SHL = self.UL - self.CL
                if S <= SHL:
                    #Sampled point is left data state
                    D = self.DL
                    U = self.UL
                    P = self.PL
                else:
                    CML = self.CL*(PM/self.PL)**self.G1;
                    STL = UM - CML
                    
                    if S > STL:
                        #Sampled point is Star Left state
                        D = self.DL*(PM/self.PL)**(1.0/self.GAMMA)
                        U = UM
                        P = PM
                    else:
                        #Sampled point is inside left fan
                        U = self.G5*(self.CL + self.G7*self.UL + S)
                        C = self.G5*(self.CL + self.G7*(self.UL - S))
                        D = self.DL*(C/self.CL)**self.G4
                        P = self.PL*(C/self.CL)**self.G3
            else:
                
                #left shock
                PML = PM/self.PL
                SL  = self.UL - self.CL*sqrt(self.G2*PML + self.G1)
                
                if S <= SL:
                    #Sampled point is left data state
                    D = self.DL
                    U = self.UL
                    P = self.PL
                    
                else:
                    #Sampled point is Star Left state
                    D = self.DL*(PML + self.G6)/(PML*self.G6 + 1.0)
                    U = UM
                    P = PM
        else:
            #Sampling point lies to the right of the contact discontinuity
            if PM > self.PR:
                #right shock
                PMR = PM/self.PR
                SR  = self.UR + self.CR*sqrt(self.G2*PMR + self.G1)
                
                if S >= SR:
                    #Sampled point is right data state
                    D = self.DR
                    U = self.UR
                    P = self.PR
                else:
                    #Sampled point is Star Right state
                    D = self.DR*(PMR + self.G6)/(PMR*self.G6 + 1.0)
                    U = UM
                    P = PM
            else:
                #Right rarefaction
                SHR = self.UR + self.CR
                if S >= SHR:
                    #Sampled point is right data state
                    D = self.DR
                    U = self.UR
                    P = self.PR
                else:
                    CMR = self.CR*(PM/self.PR)**self.G1
                    STR = UM + CMR
                    
                    if S <= STR:
                        #Sampled point is Star Right state
                        D = self.DR*(PM/self.PR)**(1.0/self.GAMMA)
                        U = UM
                        P = PM
                        
                    else:
                        #Sampled point is inside left fan
                        U = self.G5*(-self.CR + self.G7*self.UR + S)
                        C = self.G5*(self.CR - self.G7*(self.UR - S))
                        D = self.DR*(C/self.CR)**self.G4
                        P = self.PR*(C/self.CR)**self.G3
                        
        return PM, UM, S, D, U, P

if __name__ == "__main__":
    DOMLEN = 2.0;
    DIAPH1 = 1.0;
    CELLS = 1000;
    GAMMA = 5.0/3.0;
    R = 1.0;
    DL = 1.0;
    UL = 0;
    PL = 1.0;
    DR = 0.125*DL;
    UR = 0.0;
    PR = 0.1*PL;
    
    
    DSCALE = DL;
    TSCALE = 1.0;
    PSCALE = DSCALE*R*TSCALE;
    USCALE = sqrt(R*TSCALE)

    TIMEOU = 0.5
    
    R = riemannSolver(DOMLEN,DIAPH1,CELLS,GAMMA, R, DL, UL, PL, DR, UR, PR, USCALE, DSCALE, PSCALE, TSCALE)
    
    R.compute(TIMEOU, plot_this = True)