import matplotlib.pyplot as plt
import numpy as np
import sympy as sym

x = sym.Symbol("x")
h = sym.Symbol("h")

N = 2
step = 1/(2*N)

psi0 = (x-h)*(x-2*h)/(2*h*h)
psi1 = -x*(x-2*h)/(h*h)
psi20 = x*(x-h)/(2*h*h)

d_psi0 = sym.simplify(sym.diff(psi0, x))
d_psi1 = sym.simplify(sym.diff(psi1, x))
d_psi20 = sym.simplify(sym.diff(psi20, x))

psi00 = sym.simplify(d_psi0*d_psi0)
psi01 = sym.simplify(d_psi0*d_psi1)
psi02 = sym.simplify(d_psi0*d_psi20)
psi11 = sym.simplify(d_psi1*d_psi1)
psi12 = sym.simplify(d_psi1*d_psi20)
psi22_0 = sym.simplify(d_psi20*d_psi20)

A00 = sym.simplify(sym.integrate(psi00, (x, 0, 2*h)))
A01 = sym.simplify(sym.integrate(psi01, (x, 0, 2*h)))
A02 = sym.simplify(sym.integrate(psi02, (x, 0, 2*h)))
A11 = sym.simplify(sym.integrate(psi11, (x, 0, 2*h)))
A12 = sym.simplify(sym.integrate(psi12, (x, 0, 2*h)))
A22_0 = sym.simplify(sym.integrate(psi22_0, (x, 0,2*h)))


psi21 = (x-1+h)*(x-1)/(2*h*h)
psi3 = -(x-1+2*h)*(x-1)/(h*h)
#psi4 = (x-1+2*h)*(x-1+h)/(2*h*h)

d_psi21 = sym.diff(psi21, x)
d_psi3 = sym.diff(psi3, x)
#d_psi4 = sym.diff(psi4, x)

psi22_1 = sym.simplify(d_psi21*d_psi21)
psi23 = sym.simplify(d_psi21*d_psi3)
#psi24 = sym.simplify(d_psi21*d_psi4)
psi33 = sym.simplify(d_psi3*d_psi3)
#psi34 = sym.simplify(d_psi3*d_psi4)
#psi44 = sym.simplify(d_psi4*d_psi4)

A22_1 = sym.simplify(sym.integrate(psi22_1, (x, 1-2*h, 1)))
A23 = sym.simplify(sym.integrate(psi23, (x, 1-2*h, 1)))
#A24 = sym.simplify(sym.integrate(psi24, (x, 1-2*h,1)))
A33 = sym.simplify(sym.integrate(psi33, (x, 1-2*h,1)))
#A34 = sym.simplify(sym.integrate(psi34, (x, 1-2*h,1)))
#A44 = sym.simplify(sym.integrate(psi44, (x, 1-2*h,1)))


A = np.array([[A00.subs(h, step), A01.subs(h, step), A02.subs(h, step),0],
                [A01.subs(h, step), A11.subs(h, step), A12.subs(h, step),0],
             [A02.subs(h, step), A12.subs(h, step), A22_0.subs(h, step) + A22_1.subs(h, step), A23.subs(h, step)],
             [0,0, A23.subs(h, step), A33.subs(h, step)]], dtype="float")

print(A)

f = 2*x - 1

b0 = sym.simplify(f*psi0)
b1 = sym.simplify(f*psi1)
b20 = sym.simplify(f*psi20)
b21 = sym.simplify(f*psi21)
b3 = sym.simplify(f*psi3)
#b4 = sym.simplify(f*psi4)

B0 = sym.simplify(sym.integrate(b0, (x, 0, 2*h)))
B1 = sym.simplify(sym.integrate(b1, (x, 0, 2*h)))
B20 = sym.simplify(sym.integrate(b20, (x, 0, 2*h)))

B21 = sym.simplify(sym.integrate(b21, (x, 1-2*h, 1)))
B3 = sym.simplify(sym.integrate(b3, (x, 1-2*h, 1)))
#B4 = sym.simplify(sym.integrate(b4, (x, 1-2*h, 1)))



def analytical(x, C, D):
    return 0.5*x*x - 1/3 * x**3 + C*x + D - C - 1/6

C = 0
D = 0


#A = np.array([[7/(6*h), -4/(3*h), 1/(6*h), 0,0,0], [-4/(3*h), 8/(3*h), -4/(3*h),0,0,0], [1/(6*h), -4/(3*h), 7/(h*6),0,0,0], [0,0,0,7/(6*h), -4/(3*h), 1/(6*h)], [0,0,0,-4/(3*h), 8/(3*h), -4/(3*h)], [0,0,0,1/(6*h), -4/(3*h), 7/(6*h)]]).reshape(6,6)



#b = np.array([-h/3, 4*h*(2*h-1)/(3), 4*h*(h - 1)/(3), h*(1-4*h)/(3), 4*h*(1-2*h)/(3), h/3]).reshape(6,1)

B = np.array([B0.subs(h, step), B1.subs(h, step), B20.subs(h, step) + B21.subs(h, step), B3.subs(h, step)], dtype="float").reshape(4,1)
#B = np.array([B0, B1, B20, B21, B3])


X = np.linalg.solve(A,B)


h = 1/(2*N)

def u(x):
    psi0 = (x-h)*(x-2*h)/(2*h*h)
    psi1 = -x*(x-2*h)/(h*h)
    psi2 = x*(x-h)/(2*h*h)

    psiNm2 = (x-1+h)*(x-1)/(2*h*h)
    psiNm1 = -(x-1+2*h)*(x-1)/(h*h)
    #psiN = (x-1+2*h)*(x-1+h)/(2*h*h)
    if x <= 0.5:
        return X[0]*psi0 + X[1]*psi1 + X[2]*psi2
    else:
        return X[2]*psiNm2 + X[3]*psiNm1 # + X[5]*psiN


u = np.vectorize(u)

x = np.linspace(0,1,1000)
plt.plot(x, analytical(x, C, D), label="exact")
plt.plot(x, u(x), label="numerical")
plt.legend()
plt.show()
