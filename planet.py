"""
A rotating planet with simple weather.
"""
import numpy as np
from scipy import ndimage

class Ring():
    """
    Parameters needed:
    r: radius
    a_t: atmosphere thickness
    a_d_in: dissipaion of incoming radiation
    a_a_in: absorbance of incoming radiation

    """
    def __init__(self, width, girth):
        """
        Parameters at each point on the surface:
        temp, pressure, vel_x, vel_y, 

        """
        self.surface = np.zeros((width, girth, 3))
        self.findheatin(width, girth)

    def rotate(self, days=1):
        for i in range(int(days*self.surface.shape[1])):
            self.heatin = np.roll(self.heatin,-1,1)
            self.heatchange()

    def orbit(self):
        pass
    
    def heatchange(self, a=1,b=1/50000,c=3):
        f = self.surface[:,:,0]
        f += self.heatin*a                  #rad in 
        f -= f**4*b                         #rad out
        f += ndimage.laplace(f)/4           #conductance

    def findheatin(self, width, girth, rad_in=1, r=1, thickness=0.2, a_d_in=1):
        """
        Heat from the sun as a function of radius (r), angle (a) 
        and atmosphere thickness (b)
        """
        angles = np.array([(i+1)*np.pi/(girth/2) for i in range(girth//2)])
        c = r**2/(1+np.tan(angles)**2)
        y_1 = np.sqrt((r+thickness)**2-c)
        y_2 = np.sqrt(c*np.tan(angles)**2)
        y = y_1 - y_2
        b = rad_in*np.sin(angles)/(y*a_d_in)
        self.heatin = np.pad(np.repeat([b],width,0),((0,0),(0,(girth+1)//2)),'constant')
        
