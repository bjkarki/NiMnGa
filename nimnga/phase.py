import numpy as np
import sys

####################################################################################################

## Class for cubic NiMnGa
class Cubic:
    
    ## Initialization
    def __init__(self, lattice):
        # Make sure the input is in the right form
        if len(lattice) == 1:
            # Lattice = a in nm
            self.lattice = lattice
        else:
            sys.exit("Insert lattice parameter of cubic crystal: a in nm")

    ## Transformation: cubic space <-> cartesian space
    ## ca_P_cu
    def cartesian(self):
        # Obtain lattice parameters
        a,b,c,Gamma = self.lattice
        # convert degree to radians.
        Gamma = np.deg2rad(Gamma)
        # Transformation to cartesian coordinates
        # properties: x||a_m, z||c_m
        caPcu = np.matrix([[a,0.,0.],\
                        [0.,a,0.], \
                        [0.,0.,a]])
        del a,b,c,Gamma
        # return 3_P_m transformation matrix
        return caPcu
    
####################################################################################################

## Class for Non-modulated NiMnGa
class Tetragonal:
    
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
class Monoclinic:
    
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