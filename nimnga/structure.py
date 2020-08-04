import sys
import numpy as np

####################################################################################################

## An instanct of NiMnGa class
class __Structure:
    
    ## Calculate angles between planes and directions
    def getangles(self,m1,m2,pord = 'd'):
        ## Structure matrix: cartesian space <-> NMG space
        Pm = self.cartesian()
        # if planes
        if pord == 'p':
            m1 = np.squeeze(np.asarray(m1 @ np.linalg.inv(Pm)))
            m2 = np.squeeze(np.asarray(m2 @ np.linalg.inv(Pm)))
        # if directions
        else:
            m1 = np.squeeze(np.asarray(Pm @ m1))
            m2 = np.squeeze(np.asarray(Pm @ m2))
        # Numerator and denominator  
        num = np.dot(m1,m2)
        den = np.linalg.norm(m1) * np.linalg.norm(m2)
        # return angle in degrees
        return np.rad2deg(np.arccos(num/den))

    ## Calculate d-spacing of rational planes & magnitudes of directions
    def getmagnitude(self,m,pord = 'd'):
        ## Structure matrix: cartesian space <-> NMG space
        Pm = self.cartesian()
        # if planes
        if pord == 'p':
            m = np.squeeze(np.asarray(m @ np.linalg.inv(Pm)))
            mag = 1/np.linalg.norm(m)
        # if directions
        else:
            m = np.squeeze(np.asarray(Pm @ m))
            mag = np.linalg.norm(m)
        # return magnitude in nm
        return mag
    
####################################################################################################

## Class for Non-modulated NiMnGa
class Tetragonal(__Structure):
    
    ## In the non-modulated NiMnGa, c > a
    
    ## Initialization
    def __init__(self, lattice):
        # Make sure the input is in the right form
        if len(lattice) == 2:
            # Lattice = [a,c] in nm
            self.lattice = lattice
        else:
            sys.exit("Insert lattice parameter of tetragonal crystal: [a, c] (nm and degree)")
    
    ## Transformation: tetragonal space <-> cartesian space
    ## c_P_t
    def cartesian(self):
        # Obtain lattice parameters
        a,c = self.lattice
        # Transformation to cartesian coordinates
        # properties: x||a_t, z||c_t
        cPt = np.matrix([[a,0.,0.],\
                        [0.,a,0.], \
                        [0.,0.,c]])
        del a,c
        # return c_P_t transformation matrix
        return cPt

####################################################################################################

## Class for 10M NiMnGa
class Monoclinic(__Structure):
    
    ## In the monoclinic 10M NiMnGa,c is the unique axis & c < b < a
    
    ## Initialization
    def __init__(self, lattice):
        # Make sure the input is in the right form
        if len(lattice) == 4:
            # Lattice = [a,b,c,Gamma] in nm and degree
            self.lattice = lattice
        else:
            sys.exit("Insert lattice parameter of monoclinic crystal: [a, b, c, gamma] (nm and degree)")

    ## Transformation: monoclinic space <-> cartesian space
    ## c_P_m
    def cartesian(self):
        # Obtain lattice parameters
        a,b,c,Gamma = self.lattice
        # convert degree to radians.
        Gamma = np.deg2rad(Gamma)
        # Transformation to cartesian coordinates
        # properties: x||a_m, z||c_m
        cPm = np.matrix([[a,b*np.cos(Gamma),0.],\
                        [0.,b*np.sin(Gamma),0.], \
                        [0.,0.,c]])
        del a,b,c,Gamma
        # return 3_P_m transformation matrix
        return cPm
    
####################################################################################################