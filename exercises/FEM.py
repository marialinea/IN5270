import matplotlib.pyplot as plt
import numpy as np
import sympy as sym
"""
A = np.array([[2,1,0,0],[1,2,0,0], [0,0,2,1], [0,0,1,2]]).reshape(4,4) * 1/12

b = np.array([31/96, 37/96, 15/32, 47/96])


X = np.linalg.solve(A,b)


def f(x):
    return 1 + 2*x - x*x

def u(x):
    phi0 = -2*x + 1
    phi10 = 2*x
    phi11 = -2*(x-1)
    phi2 = 2*x - 1
    if x <= 0.5:
        return X[0]*phi0 + X[1]*phi10
    else:
        return  X[2]*phi11 + X[3]*phi2


u = np.vectorize(u)

x = np.linspace(0,1,1000)
plt.plot(x, f(x), label="exact")
plt.plot(x, u(x), label="numerical")
plt.legend()
plt.show()
"""

N = 2
step = 1/(2*N)
x = sym.Symbol("x")
h = sym.Symbol("h")

f = 1 + 2*x - x*x

psi0 = (x-h)*(x-2*h)/(2*h*h)
psi1 = -x*(x-2*h)/(h*h)
psi20 = x*(x-h)/(2*h*h)

psi21 = (x-1+h)*(x-1)/(2*h*h)
psi3 = -(x-1+2*h)*(x-1)/(h*h)
psi4 = (x-1+2*h)*(x-1+h)/(2*h*h)


psi00 = sym.simplify(psi0*psi0)
psi01 = sym.simplify(psi0*psi1)
psi02 = sym.simplify(psi0*psi20)
psi11 = sym.simplify(psi1*psi1)
psi12 = sym.simplify(psi1*psi20)
psi22_0 = sym.simplify(psi20*psi20)

A00 = sym.simplify(sym.integrate(psi00, (x, 0, 2*h)))
A01 = sym.simplify(sym.integrate(psi01, (x, 0, 2*h)))
A02 = sym.simplify(sym.integrate(psi02, (x, 0, 2*h)))
A11 = sym.simplify(sym.integrate(psi11, (x, 0, 2*h)))
A12 = sym.simplify(sym.integrate(psi12, (x, 0, 2*h)))
A22_0 = sym.simplify(sym.integrate(psi22_0, (x, 0,2*h)))

A1 = np.array([[A00.subs(h, step), A01.subs(h, step), A02.subs(h, step)], [A01.subs(h, step), A11.subs(h, step), A12.subs(h, step)], [A02.subs(h, step), A12.subs(h, step), A22_0.subs(h, step)]]).reshape(3,3)

A1 = np.array([[A00, A01, A02], [A01, A11, A12], [A02, A12, A22_0]]).reshape(3,3)
print(A1)

psi22_1 = sym.simplify(psi21*psi21)
psi23 = sym.simplify(psi21*psi3)
psi24 = sym.simplify(psi21*psi4)

psi33 = sym.simplify(psi3*psi3)
psi34 = sym.simplify(psi3*psi4)

psi44 = sym.simplify(psi4*psi4)

A22_1 = sym.simplify(sym.integrate(psi22_1, (x, 2*h, 1)))
A23 = sym.simplify(sym.integrate(psi23, (x, 2*h, 1)))
A24 = sym.simplify(sym.integrate(psi24, (x, 2*h,1)))

A33 = sym.simplify(sym.integrate(psi33, (x, 2*h,1)))
A34 = sym.simplify(sym.integrate(psi34, (x, 2*h,1)))
A44 = sym.simplify(sym.integrate(psi44, (x, 2*h,1)))

A2 = np.array([[A22_1.subs(h, step), A23.subs(h, step), A24.subs(h, step)], [A23.subs(h, step), A33.subs(h, step), A34.subs(h, step)], [A24.subs(h, step), A34.subs(h, step), A44.subs(h, step)]]).reshape(3,3)



#A = np.block([[A1, zero], [zero, A2]])
#A = np.array(A)
A = np.array([[A00.subs(h, step), A01.subs(h, step), A02.subs(h, step),0,0,0], [A01.subs(h, step), A11.subs(h, step), A12.subs(h, step),0,0,0], [A02.subs(h, step), A12.subs(h, step), A22_0.subs(h, step), 0,0,0], [0,0,0, A22_1.subs(h, step), A23.subs(h, step), A24.subs(h, step)], [0,0,0, A23.subs(h, step), A33.subs(h, step), A34.subs(h, step)], [0,0,0, A24.subs(h, step), A34.subs(h, step), A44.subs(h, step)]], dtype="float")


b0 = sym.simplify(f*psi0)
b1 = sym.simplify(f*psi1)
b20 = sym.simplify(f*psi20)
b21 = sym.simplify(f*psi21)
b3 = sym.simplify(f*psi3)
b4 = sym.simplify(f*psi4)

B0 = sym.simplify(sym.integrate(b0, (x, 0, 2*h)))
B1 = sym.simplify(sym.integrate(b1, (x, 0, 2*h)))
B20 = sym.simplify(sym.integrate(b20, (x, 0, 2*h)))
B21 = sym.simplify(sym.integrate(b21, (x, 2*h, 1)))
B3 = sym.simplify(sym.integrate(b3, (x, 2*h, 1)))
B4 = sym.simplify(sym.integrate(b4, (x, 2*h, 1)))

B = np.array([B0.subs(h, step), B1.subs(h, step), B20.subs(h, step), B21.subs(h, step), B3.subs(h, step), B4.subs(h, step)], dtype="float").reshape(6,1)




X = np.linalg.solve(A,B)
print(X)
def u(x,h):
    psi0 = (x-h)*(x-2*h)/(2*h*h)
    psi1 = -x*(x-2*h)/(h*h)
    psi20 = x*(x-h)/(2*h*h)

    psi21 = (x-1+h)*(x-1)/(2*h*h)
    psi3 = -(x-1+2*h)*(x-1)/(h*h)
    psi4 = (x-1+2*h)*(x-1+h)/(2*h*h)
    if x <= 0.5:
        return X[0]*psi0 + X[1]*psi1 + X[2]*psi20
    else:
        return  X[3]*psi21 + X[4]*psi3 + X[5]*psi4


def f(x):
    return 1 + 2*x - x*x


u = np.vectorize(u)

x = np.linspace(0,1,1000)
plt.plot(x, f(x), label="exact")
plt.plot(x, u(x,step), label="numerical")
plt.legend()
plt.show()
