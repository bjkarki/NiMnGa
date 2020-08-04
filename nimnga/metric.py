import numpy as np
import sys

####################################################################################################

## Calculate angles between planes and directions
def angles(self,m1,m2,pord = 'd'):
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

####################################################################################################

## Calculate d-spacing of rational planes & magnitudes of directions
def magnitude(self,m,pord = 'd'):
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