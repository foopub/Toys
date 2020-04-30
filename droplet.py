"""
Create and optimise a droplet with surface area A, volume V,
density D and potentially add a surface tension factor.
"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy import optimize


def calc_area(y: np.array) -> float:
    """
    Calculate area of revolution using Pappus's centroid theorem.
    A = sd, 
    s = sqrt(dx^2+dy^2), 
    d = 2pi*c = 2pi*(y+dy/2) = pi(2y+dy)
    """
    dx = y[-1]
    dy = y[1:-1]-y[:-2]
    s = np.sqrt(dy**2+dx**2)
    d = 2*np.pi*(y[:-2]+dy/2)
    return np.sum(s*d)

def calc_vol(y: np.array) -> np.array:
    """
    Calculate volume of revolution.

    V = Ad          #A is area, d is centroid travel distance  
    A = h*dx,                       #h is mean height
    d = 2pi*c = 2pi*h/2 = pi*h      #c is the centroid y coord.
    h = (y_0+y_1)/2
    V = h*dx*pi*h = h^2*pi*dx
    """
    dx = y[-1]
    h = (y[1:-1]+y[:-2])/2
    return np.sum(np.pi*(h**2)*dx)

def calc_pe(y: np.array, x_points: np.array) -> float:
    """
    Calculate relative potential energy.
    PE_r = vol*x
    """
    dx = y[-1]
    h = (y[1:-1]+y[:-2])/2
    vol = np.pi*(h**2)*dx
    return -np.sum(vol*x_points)

def constraint_func(y: np.array) -> np.array:
    """
    Combine the constraints into a single array.
    """
    return np.array([calc_area(y),calc_vol(y)])

def optimise(
        res: int=10,
        area: int=10,
        radius: int=1,
        volume: int=2, 
        ) -> optimize.OptimizeResult:

    assert volume > (area-np.pi*radius**2)**(3/2)/(6*np.sqrt(np.pi))
    func_ub = np.array([area,volume])
    func_lb = np.array([area,volume])
    nonlinear_constraint = optimize.NonlinearConstraint(
            constraint_func, func_lb, func_ub)
    bounds = optimize.Bounds([radius,*np.zeros(res+1)],
            [radius,*[np.inf for _ in range(res-1)],0,np.inf])

    y = np.linspace(radius,0,res+1)
    y0 = [*y,0.05]
    x_points = np.arange(0.5,res,1)

    result = optimize.minimize(calc_pe,y0,args=(x_points),
            method='trust-constr', bounds=bounds,
            constraints=[nonlinear_constraint],
            options={'verbose': 1,'maxiter': 3000})
                #'xtol': 1e-13, 'gtol': 1e-13})

    r = result.x[:-1]
    z = [i*result.x[-1] for i in range(len(result.x)-1)]
    #theta = np.r_[0:2*np.pi:30j]
    #x = r*np.sin(theta)
    plt.plot(z,r)
    return result 



