import sys
import numpy as np

####################################################################################################

## twinning mode in 10M NiMnGa

from . import structure # monoclinic phase required

class Monoclinic(structure.Monoclinic):
    
    ## Initialization
    def __init__(self,lattice):
        super().__init__(lattice)
    
    ## function to calculate irrational elements of type I and type II twins
    def __irrational(self):
        # Obtain lattice parameters
        a,b,c,Gamma = self.lattice
        # proper scaling: [Karki, Pond, et. al]
        beta = b/a
        gamma = c/a
        Gamma = np.deg2rad(Gamma)
        # irrational elements for form: k2 = (q1 -1 1), gamma1 = [r1 -1 1]
        q1 = 2*beta*np.cos(Gamma)/(gamma**2 - beta**2)
        r1 = 2*beta*gamma**2 * np.cos(Gamma) / (gamma**2 - (beta*np.sin(Gamma))**2)
        # irrational elements for form: k2 = (-1 q2 1), gamma1 = [-1 r2 1]
        q2 = q2 = 2*beta*np.cos(Gamma)/(gamma**2 - 1)
        r2 = 2*gamma**2*np.cos(Gamma)/(beta*(gamma**2 - np.sin(Gamma)**2))
        # return output variables
        return q1, r1, q2, r2
    
    ## cubic <-> monoclinic: 12 possible twinning modes 
    def twinmode(self):
        # Obtain irrational elements
        q1,r1,q2,r2 = self.__irrational()
        # 6 possiblities for k1
        k1 = ([1,1,0],[0,1,0],[0,1,1],[0,-1,1],[1,0,1],[-1,0,1])
        # 6 possibilites for k2
        k2 = ([-1,1,0],[1,0,0],[q1,-1,1],[-q1,1,1],[-1,q2,1],[1,-q2,1])
        # 6 possiblities for gamma1
        gamma1 = ([-1,1,0],[1,0,0],[r1,-1,1],[-r1,1,1],[-1,r2,1],[1,-r2,1])
        # 6 possibilities for gamma2
        gamma2 = ([1,1,0],[0,1,0],[0,1,1],[0,-1,1],[1,0,1],[-1,0,1])
        # Return twinning modes
        return k1, k2, gamma1, gamma2
    
    ## 2' about a direction or plane normal in cartesian coordinates space
    def __pirotation (self,m,pord = 'd'):
        # Obtain 3_P_m transformation matrix
        Pm = self.cartesian()
        # if planes
        if pord == 'p':
            m = m @ np.linalg.inv(Pm)
        # if directions
        else:
            m = Pm @ m
        # normalize the vector
        m = m/np.linalg.norm(m) # crucial step!
        m = np.squeeze(np.asarray(m))
        # check if the first argument has correct size
        if len(m) == 3:
            m1, m2, m3 = m
        else:
            sys.exit("The vector should be 3 dimensional")
        # 2pi rotation about the unit normal m
        R = np.asarray([[2*m1**2-1, 2*m1*m2, 2*m1*m3],\
                        [2*m2*m1, 2*m2**2-1, 2*m2*m3],\
                        [2*m3*m1, 2*m3*m2, 2*m3**2-1]])
        # return output variables
        return R
        
####################################################################################################
