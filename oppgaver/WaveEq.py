import OneDim_WaveSolver as WaveSolver
import numpy as np


Nx = [10**i for i in range(1,4)]
L = 1
c = 1
T = 1

def I(x):
    return np.sin((2*np.pi*x)/L)

def E(x,t):
    return np.sin((2*np.pi*x)/L)*np.cos((2*np.pi*c*t)/L)

error = np.zeros(len(Nx))
stepsize = np.zeros(len(Nx))

for i in range(len(Nx)):
    print("Running for Nx = {}".format(Nx[i]))
    MySolver = WaveSolver.OneDim_WaveSolver(Nx[i], L, c, T)
    MySolver.Conditions(I)
    MySolver.FirstTimeStep()
    MySolver.AdvanceTime(E)
    error[i] = MySolver.l_infty_norm
    stepsize[i] = MySolver.dt

for i in range(2):
    r = np.log(error[i+1]/error[i])/np.log(stepsize[i+1]/stepsize[i])
    print(r)

# Konvergensraten er feil
