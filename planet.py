"""
Work in progress. 
A rotating planet with simple weather. Using this to learn simulation
techniques and nonlinear dynamics, etc. 
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import ndimage

def curl(array: np.array) -> np.array:
    """
    """
    pass

def div(array: np.array) -> np.array:
    """
    Divergence of array whose first axis contains the dimensions:
    a[0] = z
    a[1] = y
    a[2] = x
    etc.
    """
    d = array.shape[0]
    return np.add.reduce([np.gradient(array[i],axis=i) for i in range(d)])

class Tube():
    def __init__(self, length: int, conductivity: float=1,
            viscosity: float=2, density: float=1.2, pressure: float=200.1):
        """
        1D fluid.
        """
        self.c = conductivity
        self.m = viscosity
        

        self.T = np.repeat(273.15, length)      #initial temperature
        self.d = np.repeat(density, length)     #densities
        self.P = np.repeat(pressure, length)    #pressure
        self.u = np.zeros(length)

        self.mass = density*length
        self.length = length

    def heatchange(self, rate_a, rate_b):
        #add some heat in
        self.T[0] += rate_a    
        #radiate some heat out
        self.T[-1] -= self.T[-1]**4*rate_b #???
        #distribute temp
        self.T -= self.c*ndimage.laplace(self.T)/2

    def velocitychange(self, dt=0.1):
        """
        Work towards full Navier Stokes equation:
        P' = P - m' div(u)
        p Du/Dt = -del(P) + pg + m lap(u) + 1/3 m del(div(u))
        
        2nd and 3rd term combine as there's only one dimension and
        del(div(u)) = lap(u)
        """
        du_dt = -np.gradient(self.P)/self.d - 9.81
        du_dt += (4/3*self.m/self.d - self.u) * ndimage.laplace(self.u)
        self.u += du_dt*dt

    def add_edges(self, array):
        array = np.insert(array,0,max(-array[0],0))
        array = np.append(array,min(-array[-1],0))
        return array

    def densitychange(self,dt=0.1):
        #Using conservation of mass
        dp = dt*-np.gradient(self.add_edges(self.u))[1:-1]
        #Change pressure accordingly to ideal gas law
        self.P *= 1+dp/self.d
        self.d += dp
        #Conserve masss by spreading out fluctuations 
        self.d += (self.mass-np.sum(self.d))/self.length

    def step(self, dt: float=0.1):
        self.velocitychange(dt)
        self.densitychange(dt)
        return np.stack([self.d, self.P, self.u])

class Sheet():
    def __init__(self, dimensions: tuple, conductivity: float=1,
            viscosity: float=2, density: float=1.2, pressure: float=200.1,
            gravity: float=9.81):
        """
        First generalisation attempt. 
        """
        self.c = conductivity
        self.m = viscosity
        self.g = gravity
        
        self.dye = np.zeros(dimensions)
        self.dye_total = 0

        self.T = np.full(dimensions, 273.15)      #initial temperature
        self.d = np.full(dimensions, density)     #densities
        self.P = np.full(dimensions, pressure)    #pressure
        #The order of coordinates is reverse:
        #For example u[:,z,y,x][0] gives the z component
        self.u = np.zeros((len(dimensions), *dimensions))

        self.mass = density*np.prod(dimensions)
        self.dim = dimensions
        self.vol = np.prod(dimensions)

        self.fig, self.ax = plt.subplots()
        self.x = np.arange(self.dim[1])
        self.y = np.arange(self.dim[0])


    def heatchange(self, rate_a: float, rate_b: float):
        #add some heat in
        self.T[0,:] += rate_a
        #radiate some heat out
        self.T -= self.T[-1]**4*rate_b #???
        #distribute temp
        self.T += ndimage.laplace(self.T)/4

    def velocitychange(self, dt=0.1):
        """
        Work towards full Navier Stokes equation:
        P' = P - m' div(u)
        p Du/Dt = -del(P) + pg + m lap(u) + 1/3 m del(div(u))
        
        """
        du_dt = -np.array(np.gradient(self.P))/self.d
    #dimensino specific du_dt[-2]
        du_dt[0] -= self.g   #gravity only in the y direction
        du_dt += self.m*np.stack(
                [ndimage.laplace(self.u[i]) for i in [0,1]])
        du_dt += 1/3*self.m*np.array(np.gradient(div(self.u)))
    #dimension specific axis=()
        #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        du_dt -= np.add.reduce(np.stack((self.u,)*2,axis=1)*np.array(
            np.gradient(self.u,axis=(1,2))))
        self.u += du_dt*dt

    def edge_velocity(self):
    #dimension specific
        #reflext x values at x edges
        self.u[1,:,0] = -self.u[1,:,1]
        self.u[1,:,-1] = -self.u[1,:,-2]
        #mirror x values at y edges 
        self.u[1,0,:] = self.u[1,1,:]
        self.u[1,-1,:] = self.u[1,-2,:]
        #mirror y values at x edges
        self.u[0,:,0] = self.u[0,:,1]
        self.u[0,:,-1] = self.u[0,:,-2]
        #mirror y values at y edges 
        self.u[0,0,:] = -self.u[0,1,:]
        self.u[0,-1,:] = -self.u[0,-2,:]

    def edge_pressure(self):
        self.P[:,0] = self.P[:,1]
        self.P[:,-1] = self.P[:,-2]
        self.P[0,:] = self.P[1,:]
        self.P[-1,:] = self.P[-2,:]

    def densitychange(self,dt=0.1):
        #Using conservation of mass
        dp = dt*-div(self.u)
        dp[[0,-1]] = dp[:,[0,-1]] = 0
        #Change pressure accordingly to ideal gas law
        self.P *= 1+dp/self.d
        self.d += dp
        #Conserve masss by spreading out fluctuations 
        self.d += (self.mass-np.sum(self.d))/self.vol

    def spread_dye(self, dt=0.1):
        dp = dt*-div(self.u*self.dye)
        self.dye += dp
        self.dye *= self.dye_total/np.sum(self.dye)
    
    def add_dye(self,x: int, y:int, spread: int, amount: float=1):
        try:
            self.dye[x-spread:x+spread,y-spread:y+spread] += amount
            self.dye_total += spread**2*amount
        except:
            print("Something's wrong, check the values!")

    def step(self, dt: float=0.1):
        self.velocitychange(dt)
        self.edge_velocity()
        self.densitychange(dt)
        self.spread_dye(dt)
        self.edge_pressure()

    def show(self, which: str="dye", arrows: bool=False):
        if which == "dye":
            p = self.ax.pcolormesh(self.dye)
        elif which == "b":
            p = self.ax.pcolormesh(self.dye)
            p = self.ax.pcolormesh(self.P)
        elif which == "p":
            p = self.ax.pcolormesh(self.P)
        else:
            p = self.ax.pcolormesh(self.d)
        if arrows:
            q = self.ax.quiver(self.x+0.5,self.y+0.5,
                    self.u[1]*self.d,self.u[0]*self.d)

    def animation(self, step: int, dt: float=0.01):
        for _ in range(10):
            self.step(dt)
        plt.cla()
        self.ax.pcolormesh(self.d)
        self.ax.quiver(self.x+0.5,self.y+0.5,
            self.u[1]*self.d,self.u[0]*self.d)

if __name__ == "__main__":
    a = Sheet((100,100),gravity=0)
    for i in range(10):
        a.u[0,i+45,abs(5-i)] += 20
    #for _ in range(100):
    #    a.step(0.01)
    ani = FuncAnimation(a.fig, a.animation, frames=800,
            interval=30, repeat=False)
    ani.save("plop4.mp4")
