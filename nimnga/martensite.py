## Import required packages
import numpy as np
import sys



## Class for 10M NiMnGa
class NMG_10M:
  
  ## Module1: Initialization
  def __init__(self, lattice):
    # Make sure the input is in the right form
    if len(lattice) == 4:
      # Lattice = [a,b,c,Gamma] in nm and degree
      self.lattice = lattice
    else:
      sys.exit("Insert lattice parameter of Monoclinic crystal in the form [a, b, c, gamma] (nm and degree)")


  ## Module2: Transformation from monoclinic space to 3-space
  def t3space(self):
    # Obtain lattice parameters
    a,b,c,Gamma = self.lattice
    # scale, such that a = 1.
    beta = b/a
    gamma = c/a
    Gamma = np.deg2rad(Gamma)
    # Transformation to 3-space aka cartesian coordinates
    # x || a_m & |x| = |a_m| = 1
    # z || c_m
    Pm = np.matrix([[1.,beta*np.cos(Gamma),0.], \
                    [0.,beta*np.sin(Gamma),0.], \
                    [0.,0.,gamma]])
    # return 3_P_m transformation matrix
    return Pm


  ## Module3: Calculate angles between planes and directions
  def angles(self,m1,m2,pord = 'd'):
    # Obtain 3_P_m transformation matrix
    Pm = self.t3space()
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


  ## Module4: Calculate d-spacing of rational planes & magnitudes of directions
  def magnitude(self,m,pord = 'd'):
    # Obtain 3_P_m transformation matrix
    Pm = self.t3space()
    # if planes
    if pord == 'p':
      m = np.squeeze(np.asarray(m @ np.linalg.inv(Pm)))
      mag = self.lattice[0]/np.linalg.norm(m)
    # if directions
    else:
      m = np.squeeze(np.asarray(Pm @ m))
      mag = np.linalg.norm(m)*self.lattice[0]
    # return magnitude in nm
    return mag


