from real_shit import P2_FEM


def f(x):
    return 2*x - 1

C = 0
D = 0
Ne = 3
symbolic = False

my_solver = P2_FEM(Ne, f, C, D, symbolic)
my_solver.BasisFunctions()
