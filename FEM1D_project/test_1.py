import matplotlib.pyplot as plt
import numpy as np
import sympy as sym

N = 1
C = 0
D = 0

x = sym.Symbol("x")
h = sym.Symbol("h")

psi0 = (x-h)*(x-2*h)/(2*h*h)
psi1 = -x*(x-2*h)/(h*h)
#psi2 = x*(x-h)/(2*h*h)

d_psi0 = sym.simplify(sym.diff(psi0, x))
d_psi1 = sym.simplify(sym.diff(psi1, x))
#d_psi2 = sym.simplify(sym.diff(psi2, x))

psi00 = sym.simplify(d_psi0*d_psi0)
psi01 = sym.simplify(d_psi0*d_psi1)
#psi02 = sym.simplify(d_psi0*d_psi2)
psi11 = sym.simplify(d_psi1*d_psi1)
#psi12 = sym.simplify(d_psi1*d_psi2)
#psi22 = sym.simplify(d_psi2*d_psi2)

A00 = sym.simplify(sym.integrate(psi00, (x, 0, 2*h)))
A01 = sym.simplify(sym.integrate(psi01, (x, 0, 2*h)))
#A02 = sym.simplify(sym.integrate(psi02, (x, 0, 1)))
A11 = sym.simplify(sym.integrate(psi11, (x, 0, 2*h)))
#A12 = sym.simplify(sym.integrate(psi12, (x, 0, 1)))
#A22 = sym.simplify(sym.integrate(psi22, (x, 0, 1)))

step = 1/(2*N)

#A1 = np.array([[A00.subs(h, step), A01.subs(h, step), A02.subs(h, step)], [A01.subs(h, step), A11.subs(h, step), A12.subs(h, step)], [A02.subs(h, step), A12.subs(h, step), A22.subs(h, step)]], dtype="float").reshape(3,3)

#A1 = np.array([[A00, A01, A02], [A01, A11, A12], [A02, A12, A22]]).reshape(3,3)

A2 = np.array([[A00.subs(h, step), A01.subs(h, step)], [A01.subs(h, step), A11.subs(h, step)]], dtype="float").reshape(2,2)

f = 2*x - 1

b0 = sym.simplify(f*psi0)
b1 = sym.simplify(f*psi1)
#b2 = sym.simplify(f*psi2)


B0 = sym.simplify(sym.integrate(b0, (x, 0, 2*h)))
B1 = sym.simplify(sym.integrate(b1, (x, 0, 2*h)))
#B2 = sym.simplify(sym.integrate(b2, (x, 0, 1)))

#B = np.array([B0.subs(h, step), B1.subs(h, step), B2.subs(h, step)], dtype="float").reshape(3,1)
B2 = np.array([B0.subs(h, step), B1.subs(h, step)], dtype="float").reshape(2,1)
#B = np.array([B0, B1, B2])
#print(B)


h = 1/(2*N)


X = np.linalg.solve(A2,B2)


def analytical(x, C, D):
    return 0.5*x*x - 1/3 * x**3 + C*x +D -C -1/6

def u(x):
    psi0 = (x-h)*(x-2*h)/(2*h*h)
    psi1 = -x*(x-2*h)/(h*h)
    #psi2 = x*(x-h)/(2*h*h)
    return X[0]*psi0 + X[1]*psi1 #+ X[2]*psi2


u = np.vectorize(u)


x = np.linspace(0,1,1000)
plt.plot(x, analytical(x, C, D), label="exact")
plt.plot(x, u(x), label="numerical")
plt.legend()
plt.show()

"""
psi0 = lambda x, h: (x-h)*(x-2*h)/(2*h*h)
psi1 = lambda x, h: -x*(x-2*h)/(h*h)
psi2 = lambda x, h: x*(x-h)/(2*h*h)




x = np.linspace(0,1,1000)


plt.plot(x, psi0(x,h), label="psi0")
plt.plot(x, psi1(x,h), label="psi1")
plt.plot(x, psi2(x,h), label="psi2")
plt.legend()
plt.show()
"""
