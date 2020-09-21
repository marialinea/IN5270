import numpy as np
from WaveSolver2D import WaveSolver2D
import matplotlib.pyplot as plt

b = 0
Nx = 50
Ny = 50
Lx = 1
Ly = 1
T = 10.


def I(x,y, I0=1, Ia=1, Im=0.5, Is=0.1):
    return I0 + Ia*np.exp(-((x-Im)/Is)**2)

def V(x,y):
    return 0

def q(x,y, H0 = 10):
    g = 9.81
    return g*(H0 - B(x,y))

"""
def B(x, y, B0=0, Ba=1, Bmx=0.5, Bmy=0.5, Bs=0.1, b=1):
    return B0 + Ba*np.exp(-((x-Bmx)/Bs)**2-((y-Bmy)/(b*Bs))**2)
"""
def B(x, y, B0=0, Ba=1, Bmx=0.5, Bmy=0.5, Bs=2, b=1):
    return B0 + Ba*np.cos(np.pi*(x-Bmx/(2*Bs)))*np.cos(np.pi*(y-Bmy/(2*Bs)))



def f(x,y,t):
    return 0



solver = WaveSolver2D(b, Nx, Ny, Lx, Ly, T)
solver.Initialize(I, V, q, f)
solver.FirstTimeStep()
solver.AdvanceTime()


X, Y = np.meshgrid(solver.x, solver.y)
plt.contourf(X.T,Y.T, solver.u[1:-1,1:-1])
plt.title("Numerical solution")
plt.xlabel("X")
plt.ylabel("Y")
plt.colorbar(label="f(x,y,$t_n$)")
plt.show()
