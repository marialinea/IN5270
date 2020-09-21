import numpy as np
from WaveSolver2D import WaveSolver2D

b = 0
Lx = 1
Ly = 1
T = 1
mx = 1; my = 1;
kx = mx*np.pi/Lx
ky = my*np.pi/Ly

def I(x,y):
    A = 1
    kx = np.pi
    ky = np.pi
    return A*np.cos(kx*x)*np.cos(ky*y)

def V(x,y):
    return 0

def q(x,y):
    return 1

def f(x,y,t):
    return 0

def analytical(x, y, t):
    A = 1
    kx = np.pi
    ky = np.pi
    w = np.sqrt(2*np.pi**2)
    return A*np.cos(kx*x)*np.cos(ky*y)*np.cos(t*w)

N = [2**i for i in range(1,8)]
h = [1/(i-1) for i in N]
error = []

for i in range(len(N)):
    solver = WaveSolver2D(b, N[i], N[i], Lx, Ly, T)
    solver.Initialize(I, V, q, f)
    solver.FirstTimeStep()
    solver.AdvanceTime()

    error.append(solver.ComputeError(analytical))

for i in range(len(N)-1):
    r = np.log10(error[i+1]/error[i])/np.log10(h[i+1]/h[i])
    print(r)
