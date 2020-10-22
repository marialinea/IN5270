import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

x = sym.Symbol("x")
h = sym.Symbol("h")
j = sym.Symbol("j")

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

#A1 = np.array([[A00.subs(h, step), A01.subs(h, step), A02.subs(h, step)], [A01.subs(h, step), A11.subs(h, step), A12.subs(h, step)], [A02.subs(h, step), A12.subs(h, step), A22_0.subs(h, step)]]).reshape(3,3)

#A1 = np.array([[A00, A01, A02], [A01, A11, A12], [A02, A12, A22_0]]).reshape(3,3)
#print(A1)

psi21 = (x-1+h)*(x-1)/(2*h*h)
psi3 = -(x-1+2*h)*(x-1)/(h*h)
psi4 = (x-1+2*h)*(x-1+h)/(2*h*h)

d_psi21 = sym.diff(psi21, x)
d_psi3 = sym.diff(psi3, x)
d_psi4 = sym.diff(psi4, x)

psi22_1 = sym.simplify(d_psi21*d_psi21)
psi23 = sym.simplify(d_psi21*d_psi3)
psi24 = sym.simplify(d_psi21*d_psi4)
psi33 = sym.simplify(d_psi3*d_psi3)
psi34 = sym.simplify(d_psi3*d_psi4)
psi44 = sym.simplify(d_psi4*d_psi4)

A22_1 = sym.simplify(sym.integrate(psi22_1, (x, 1-2*h, 1)))
A23 = sym.simplify(sym.integrate(psi23, (x, 1-2*h, 1)))
A24 = sym.simplify(sym.integrate(psi24, (x, 1-2*h,1)))
A33 = sym.simplify(sym.integrate(psi33, (x, 1-2*h,1)))
A34 = sym.simplify(sym.integrate(psi34, (x, 1-2*h,1)))
A44 = sym.simplify(sym.integrate(psi44, (x, 1-2*h,1)))

#A2 = np.array([[A22_1.subs(h, step), A23.subs(h, step), A24.subs(h, step)], [A23.subs(h, step), A33.subs(h, step), A34.subs(h, step)], [A24.subs(h, step), A34.subs(h, step), A44.subs(h, step)]]).reshape(3,3)

#A2 = np.array([[A22_1, A23, A24], [A23, A33, A34], [A24, A34, A44]]).reshape(3,3)
#print(A2)

A = np.array([[A00.subs(h, step), A01.subs(h, step), A02.subs(h, step),0,0,0], [A01.subs(h, step), A11.subs(h, step), A12.subs(h, step),0,0,0], [A02.subs(h, step), A12.subs(h, step), A22_0.subs(h, step), 0,0,0], [0,0,0, A22_1.subs(h, step), A23.subs(h, step), A24.subs(h, step)], [0,0,0, A23.subs(h, step), A33.subs(h, step), A34.subs(h, step)], [0,0,0, A24.subs(h, step), A34.subs(h, step), A44.subs(h, step)]], dtype="float")


print(A)

psiJm1 = (x-j*h)*(x-j*h-h)/(2*h*h)
psiJ= -(x-j*h+h)*(x-j*h-h)/(h*h)
psiJp1 = (x-j*h)*(x-j*h+h)/(2*h*h)

d_psiJm1 = sym.simplify(sym.diff(psiJm1, x))
d_psiJ = sym.simplify(sym.diff(psiJ, x))
d_psiJp1 = sym.simplify(sym.diff(psiJp1, x))

j_psi00 = sym.simplify(d_psiJm1*d_psiJm1)
j_psi01 = sym.simplify(d_psiJm1*d_psiJ)
j_psi02 = sym.simplify(d_psiJm1*d_psiJp1)
j_psi11 = sym.simplify(d_psiJ*d_psiJ)
j_psi12 = sym.simplify(d_psiJ*d_psiJp1)
j_psi22_0 = sym.simplify(d_psiJp1*d_psiJp1)

j_A00 = sym.simplify(sym.integrate(j_psi00, (x, (j-1)*h, (j+1)*h)))
j_A01 = sym.simplify(sym.integrate(j_psi01, (x, (j-1)*h, (j+1)*h)))
j_A02 = sym.simplify(sym.integrate(j_psi02, (x, (j-1)*h, (j+1)*h)))
j_A11 = sym.simplify(sym.integrate(j_psi11, (x, (j-1)*h, (j+1)*h)))
j_A12 = sym.simplify(sym.integrate(j_psi12, (x, (j-1)*h, (j+1)*h)))
j_A22_0 = sym.simplify(sym.integrate(j_psi22_0, (x, (j-1)*h,(j+1)*h)))

j_A1 = np.array([[j_A00, j_A01, j_A02], [j_A01, j_A11, j_A12], [j_A02, j_A12, j_A22_0]]).reshape(3,3)
#print(j_A1)



f = 2*x - 1

b0 = sym.simplify(f*psi0)
b1 = sym.simplify(f*psi1)
b20 = sym.simplify(f*psi20)
b21 = sym.simplify(f*psi21)
b3 = sym.simplify(f*psi3)
b4 = sym.simplify(f*psi4)

B0 = sym.simplify(sym.integrate(b0, (x, 0, 2*h)))
B1 = sym.simplify(sym.integrate(b1, (x, 0, 2*h)))
B20 = sym.simplify(sym.integrate(b20, (x, 0, 2*h)))
B21 = sym.simplify(sym.integrate(b21, (x, 1-2*h, 1)))
B3 = sym.simplify(sym.integrate(b3, (x, 1-2*h, 1)))
B4 = sym.simplify(sym.integrate(b4, (x, 1-2*h, 1)))

print(B20)
"""
B = np.array([B0.subs(h, step), B1.subs(h, step), B20.subs(h, step), B21.subs(h, step), B3.subs(h, step), B4.subs(h, step)], dtype="float").reshape(6,1)
#B = np.array([B0, B1, B20, B21, B3, B4])
#print(B)

X = np.linalg.solve(A,B)

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


def analytical(x, C, D):
    return 0.5*x*x - 1/3 * x**3 + C*x + D - C - 1/6


u = np.vectorize(u)

x = np.linspace(0,1,1000)
plt.plot(x, analytical(x, 0,0), label="exact")
plt.plot(x, u(x,step), label="numerical")
plt.legend()
plt.show()



"""
"""
psi0 = lambda x, h: (x-h)*(x-2*h)/(2*h*h)
psi1 = lambda x, h: -x*(x-2*h)/(h*h)
psi2 = lambda x, h: x*(x-h)/(2*h*h)

psiNm2 = lambda x, h: (x-1+h)*(x-1)/(2*h*h)
psiNm1 = lambda x, h: -(x-1+2*h)*(x-1)/(h*h)
psiN = lambda x, h: (x-1+2*h)*(x-1+h)/(2*h*h)

h = 1/4
x = np.linspace(0,1/2,1000)
y = np.linspace(0.5,1,1000)

plt.plot(x, psi0(x,h), label="psi0")
plt.plot(x, psi1(x,h), label="psi1")
plt.plot(x, psi2(x,h), label="psi2")
plt.axvline(x = 0.5)

plt.plot(y, psiNm2(y,h), label="psiNm2")
plt.plot(y, psiNm1(y,h), label="psiNm1")
plt.plot(y, psiN(y,h), label="psiN")

plt.legend()
plt.show()
"""
