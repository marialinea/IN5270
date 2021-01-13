from FEM_class import P2_FEM
import numpy as np


Ne = 3
f = "2*x-1"
analytical = "0.5*x*x -1/3*x**3 + C*x - C + D - 1/6"

N = 3
C = [10**i for i in range(-1,2)]
D = [10**i for i in range(-1,2)]

N = 31
Ne = np.arange(1,N)
h = []
L2 = []
"""
for i in Ne:
    print("i = ", i)
    my_solver = P2_FEM(i, f, 3, 2, analytical)
    my_solver.SolveLinearSystem()
    my_solver.L2Norm()

    h.append(1/i)
    L2.append(my_solver.L2)


r = []
for i in range(N-2):
    convergence = np.log10(L2[i+1]/L2[i])/np.log10(h[i+1]/h[i])
    r.append(convergence)

h_L2 = np.array([h, L2])

np.save("L2_C_3_D_2.npy", h_L2)
"""

my_solver = P2_FEM(3, f, 0, 0, analytical)
my_solver.BasisFunctions()
my_solver.MatrixAssembly()
