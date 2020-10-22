import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

def basis(d, point_distribution="uniform", symbolic=True):
    """
    Return all local basis function phi as functions of the local point X in a 1D element with d+1 nodes.
    If symbolic=True, return symbolic expressions, else return Python functions of X.
    point_distribution can be ’uniform’ or ’Chebyshev’. """
    X = sym.symbols("X")
    if d == 0:
        phi_sym = [1]
    else:
        if point_distribution == "uniform":
            if symbolic:
                # Compute symbolic nodes
                h = sym.Rational(1, d) # node spacing
                nodes = [i*h for i in range(d+1)]
            else:
                nodes = np.linspace(0, 1, d+1)
        elif point_distribution == "Chebyshev": # Just numeric nodes
            nodes = Chebyshev_nodes(0, 1, d)
        phi_sym = [Lagrange_polynomial(X, r, nodes) for r in range(d+1)]
    # Transform to Python functions
    phi_num = [sym.lambdify([X], phi_sym[r], modules="numpy") for r in range(d+1)]

    return phi_sym if symbolic else phi_num

def Lagrange_polynomial(x, i, points):
    p=1
    for k in range(len(points)):
        if k != i:
            p *= (x - points[k])/(points[i] - points[k])
    return p


def element_matrix(phi, Omega_e, symbolic=True):
    n = len(phi)
    A_e = sym.zeros(n, n)
    X = sym.Symbol("X")
    if symbolic:
        h = sym.Symbol("h")
    else:
        h = Omega_e[1] - Omega_e[0]
    detJ = 1 # dx/dX
    for r in range(n):
        for s in range(r, n):
            A_e[r,s] = sym.integrate(phi[r]*phi[s]*detJ, (X, 0, 1))
            A_e[s,r] = A_e[r,s]
    return A_e


def element_vector(f, phi, Omega_e, symbolic=True):
    n = len(phi)
    b_e = sym.zeros(n, 1)
    # Make f a function of X
    X = sym.Symbol("X")
    if symbolic:
        h = sym.Symbol("h")
    else:
        h = Omega_e[1] - Omega_e[0]
    x = (Omega_e[0] + Omega_e[1])/2 + h/2*X # mapping
    f = f.subs("x", x) # substitute mapping formula for x
    detJ = 1 # dx/dX
    for r in range(n):
        b_e[r] = sym.integrate(f*phi[r]*detJ, (X, 0, 1))
    return b_e

x = sym.Symbol("x")
X = sym.Symbol("X")
f = 1 + 2*x - x*x

phi = basis(2)
#d_phi = [sym.diff(X, phi[r]) for r in range(len(phi))]

print(phi)



A = element_matrix(phi, Omega_e=[0, 0.5, 1], symbolic=False)
b = element_vector(f,phi, Omega_e=[0, 0.5, 1], symbolic=False)


A = np.array(A).astype(np.float64)
b = np.array(b).astype(np.float64)

print(A)
print(b)
y = np.linalg.solve(A,b)


phi0 = sym.lambdify(X, phi[0], "numpy")
phi1 = sym.lambdify(X, phi[1], "numpy")
phi2 = sym.lambdify(X, phi[2], "numpy")


def u(x):
    return y[0]*phi0(x) + y[1]*phi1(x) + y[2]*phi2(x)

#def f(x, C, D):
    #return 0.5*x*x - 1/3 * x**3 + C*x + D - C - 1/6

def f(x, C, D):
    return 1 + 2*x - x*x

u = np.vectorize(u)

x = np.linspace(0,1,1000)

plt.plot(x, f(x,0,0), label="exact")
plt.plot(x, u(x), label="numerical")
plt.legend()
plt.show()
