# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 15:29:18 2020

@author: Jeremy Raskop

Based on https://doi.org/10.1364/OE.22.028427
"""

import numpy as np
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fmin 
import pandas as pd    


"""
TODO : 
    - Check that n is not n-1 by mistake
    - Beware of indices
    - To go faster, maybe there is no need to calculate the arctan (small angles)
    - Index of refraction: see that it is correct in matlab for Omega
    - Modify the critical angle to be always good for r_k<= 4 um, so that only L_k = L_lower counts
    - Make sure initial simplex is big enough (make problem adimentional?)
    - Check whether there is a mistake that gives a factor of 2 on the taper angle (or critical angle)
"""



class nf:
    def __init__(self, target_r):
        """
        Parameters
        ----------
        target_r : float
            Target radius in um
        
        self.n : contains the total number of iterations required to attain the 
            target radius

        """
        self.optimized_Lk = []
        self.optimized_nf_length = 0
        self.iterN = 0
        self.target_r = target_r
        # self.v1 = 15 # um/s
        # self.v2 = 750 # um/s
        self.v1 = 20 
        self.v2 = 1000
        self.radius_max = 62.5
        self.n, self.radii = self.get_nk(target_r, printout=False) 
        
        # The next lines import the critical angle for the adiabaticity criteria
        self.omegas_crit = np.loadtxt('critical_angle_780.txt')[30:]
        radii_crit = np.loadtxt('diameters.txt')/2
        self.radii_crit = radii_crit[30:]
        self.f = interp1d(self.radii_crit, self.omegas_crit)
        self.crit_angles=[]
        for radius in self.radii:
            if radius>62.5 or radius<0.50962:
            # print('Radius out ouf bounds at '+str(radius)+ ': Radius_max = 62.5 um, Radius_min = 0.50962 um')
                self.crit_angles.append(1.0)
            else:
                self.crit_angles.append(self.f(radius))

    def radius_km1 (self, rk):
        v1 = self.v1
        v2 = self.v2 
        return np.sqrt((v2+v1)/(v2-v1)) * rk
    
    def radius_kp1 (self, rk):
        v1 = self.v1
        v2 = self.v2 
        return np.sqrt((v2-v1)/(v2+v1)) * rk
    
    def get_nk (self, target_radius, printout = False):
        """
        Gets the number of iterations needed to reach a nanofiber with a radius target_radius
        
        The single mode fiber is considered to have a radius of 62.5 um
        
        This function starts with the last iteration target_radius, and goes back to single mode fiber
    
        Parameters
        ----------
        target_radius : target radius
            Corresponds to target radius.
    
        Returns
        -------
        Total number of iterations
    
        """
        k_total = 0 # contains the number of iterations, n
        r_momomode = 62.5 # radius of the monomode fiber in um
        radius = target_radius
        radii=[]
        while radius < r_momomode : 
            k_total += 1 
            radius = self.radius_km1(radius)
            radii.insert(0,radius)
            if printout:
                print (radius)
        # print("n = "+str(k_total))
        return k_total-2, radii[2:]
    
    def thetas (self,lag): 
        """
        Local taper angle at k'th iteration
        
        Returns a vector containing self.n values, evaluated using the lagrange polynome and the radii

        """
        tanthetas = [0,0,0]
        for k in range(3,self.n):
            tantheta_k = (self.radii[k-3] - self.radii[k-1])/((self.v2+self.v1)/self.v2* lag[k-1] + (self.v1/self.v2-1)*lag[k])
            tanthetas.append(tantheta_k)
        return np.arctan(tanthetas)  # maybe no need for arctan with such small angles
    
    def isadiabatic (self, lag, F):
        '''
        Adiabaticity criteria (2). Returns True if is adiabatic, False otherwise

        '''
        thetas = self.thetas(lag)

        for i,theta in enumerate(thetas):
            # print(i)
            if theta> F * self.crit_angles[i]:
                print(str(theta)+' > '+str(F)+'*'+str( self.crit_angles[i] ))
                return False, thetas.sum()
        return True, 0
            
        
        # TODO complete here
        
    
    def isendofheated(self,lag):
        """ 
        For a given lagrange polynomial, finds wether criteria (6) is valid for all k
        Return True if criteria is observed, False otherwise
        """
        kmax = self.n
        v1 = self.v1
        v2 = self.v2
        for k in range(kmax-1):
            if lag[k+1]>=(v2+v1)/(v2-v1) * lag[k]:
                return False
        return True
            
    
    def lagrange(self,K, L):
        """
        From the nine values L_k contained in L, returns the lagrange interpolating
        polinomial that goes through them.

        Parameters
        ----------
        L : vector
            Contains the L_k.
        K : vector
            Contains the k 

        Returns
        -------
        A function that is the Lagrange polynome.

        """
        return lagrange(K, L)
    
    def lagrange_Lmin (self, K, L):
        x = np.arange(0,self.n,1)
        lag = lagrange(K, L)
        Lk = lag(x)
        
        radii = np.array(self.radii)

        for k, r in enumerate(radii):
            if r<4:
                Lk[k] = 400
        
        # print(Lk)
        # if self.i_plot%300 == 0:
        #     self.i_plot+=1
        #     plt.title('Iteration N '+str(self.iterN))
        #     plt.clf()
        #     plt.plot(x,Lk)
        #     # plt.draw()
        #     # plt.show()
        
        return Lk  
           
        
    
    def Ltotal (self, K, L, lag):
        v1 = self.v1
        v2 = self.v2
        summ = 0
        # print(lag)
        for i in range(self.n):
            summ += lag[i]
        res =  (v2+v1)/v2 * L[0] + 2*v1/v2 * summ 
        print("Evaluated function: "+str(res))
        return res
        
    
    # def critical_angle_scalar(self,radius): # Radius in micrometers
    #     omegas = np.loadtxt('./NFiberPulling/critical_angle_780.txt')[30:]
    #     radii = np.loadtxt('./NFiberPulling/diameters.txt')/2
    #     radii = radii[30:]
    #     f = interp1d(radii, omegas)
    #     if radius>62.5 or radius<0.50962:
    #         # print('Radius out ouf bounds at '+str(radius)+ ': Radius_max = 62.5 um, Radius_min = 0.50962 um')
    #         return 1.0
    #     return f(radius)
    
    # def critical_angle(self, radius): # vectorization of critical_angle_scalar
    #     vfunc = np.vectorize(self.critical_angle_scalar)
    #     return vfunc(radius)
    
            
    def function_to_minimize(self, L, K, F, printout=False):
        '''
        This is the function that should be provided to the mimimizer
        
        The contraints are taken into account by returning very large values when they are not obeyed

        Parameters
        ----------
        K : vector
            Chosen k's though which we choose to optimize the function. Specifically, we consider a 
            Lagrange interpolating polynomial, which passes through nine points that are fixed with even 
            intervals on the scan-number axis k and are variable on the L axis, as L_k.
        L : vector
            Length of the k'th brush length contained in K
        F : float
            adiabaticity factor

        Returns
        -------
        float
            Total length of the tapered fiber if constraints are obeyed
            
            Larger values otherwise

        '''
        # lag = self.lagrange(K, L)
        self.iterN += 1
        print('Iteration N '+str(self.iterN))
        lag = self.lagrange_Lmin(K, L)
        for el in L:
            if el >= 40000 :
                if printout:
                    print(str(el)+" is too big")
                return 10**8*el
            elif el<=400:
                if printout:
                    print(str(el)+" is too small")
                return 10**8*np.abs(el)
        
        if self.isendofheated(lag):
            if printout:
                print("Torch passes the end of the heated region")
            return 10**8
        isadiab, angle = self.isadiabatic(lag, F)
        if not isadiab:
            if printout:
                print("Is not adiabatic")
            return 10**8*np.abs(angle)  # TODO: maybe put the difference between the angle and the critical angle instead
        print(L) 
        self.optimized_Lk = L
        ltot = self.Ltotal(K, L, lag)
        self.optimized_nf_length = ltot
        return ltot
    
    def mimimize(self):
        F = 0.20 
        data = np.array(pd.read_excel("Lk.xlsx"))
        K0 = np.array([1,21,41,61,81,101,121,141,161])
        Lk0 = list(np.transpose(data[K0-1])[1])
        # Lk0 = [6000,8500,8000,7500,3500,1500,1300,1100,1000]
        Lk0 = [6320, 9830, 10065, 6452, 2919, 1280, 834, 487, 350]
        self.optimized_Lk = Lk0 
        #Lk0 = [6000,6500,6700,7500,3500,1500,1300,1100,1000]
        K0 = np.array([1,21,41,61,81,101,121,141,161])
        
        ls = [self.function_to_minimize(K0,Lk0,F)]
        
        x = np.arange(0,self.n,1)
        plt.plot(x,self.lagrange_Lmin(K0,Lk0))
        
        
        fmin(self.function_to_minimize, Lk0, args=(K0,F,),
            xtol=0.0001, ftol=1, maxiter=None, maxfun=None, full_output=0, 
            disp=1, retall=0, callback=None, initial_simplex=None)
        plt.plot(x,self.lagrange_Lmin(K0,self.optimized_Lk), label='%.2f mm'%(self.optimized_nf_length/1000))
        ls.append(self.optimized_nf_length)
        
        for i in range(10):
            
            fmin(self.function_to_minimize, self.optimized_Lk, args=(K0,F,),
                 xtol=0.0001, ftol=1, maxiter=None, maxfun=None, full_output=0, 
                 disp=1, retall=0, callback=None, initial_simplex=None)
            plt.plot(x,self.lagrange_Lmin(K0,self.optimized_Lk), label='%.2f mm'%(self.optimized_nf_length/1000))
            ls.append(self.optimized_nf_length)
        print(ls)
        plt.legend()
        plt.xlim([-1, 140])
        return self.optimized_Lk
    
    
def main ():
    nano = nf(0.3)
    
    a = nano.mimimize()
    
    
    
    
    # a = nano.critical_angle(20)
    
    # print(a)
    # x = np.linspace(0.51, 60, 1000)
    # plt.plot(x,nano.critical_angle(x))

    # L = [6000,8500,8000,7500,3500,1500,1300,1100,1000]
    # L = [4000 for i in range(9)]
    # K = [1,21,41,61,81,101,121,141,161]
    # nano.lagrange_Lmin(K, L)

    # a = nano.function_to_minimize(datax, datay, 0.2)
    # print(a)
    # print(nano.Ltotal(datax, datay))
    
    # lg = nano.lagrange(datax,datay)
    # x = np.arange(0,267,1)

    # #plt.plot(x, lg(x))
    # #plt.plot(datax,datay,'o')
    
    # a = nano.thetas(lg)
    # print(a)
    # plt.plot(x,a)
    
    return a
    
if __name__ == "__main__":
    a = main()
    
    
    
    
    
    
    
    
    
    
    