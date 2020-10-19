import sympy as sym
import numpy as np

x = sym.Symbol("x")
h = sym.Symbol("h")
C = sym.Symbol("C")
D = sym.Symbol("D")

f = 2*x - 1

#psi0 = (x-h)*(x-2*h)/(2*h*h)
#psi1 = -x*(x-2*h)/(h*h)
#psi2 = x*(x-h)/(2*h*h)

psi0 = (x-1+h)*(x-1)/(2*h*h)
psi1 = -(x-1+2*h)*(x-1)/(h*h)
psi2 = (x-1+2*h)*(x-1+h)/(2*h*h)

d_psi0 = sym.diff(psi0, x)
d_psi1 = sym.diff(psi1, x)
d_psi2 = sym.diff(psi2, x)

d_psi00 = sym.simplify(d_psi0*d_psi0)
d_psi01 = sym.simplify(d_psi0*d_psi1)
d_psi02 = sym.simplify(d_psi0*d_psi2)

d_psi11 = sym.simplify(d_psi1*d_psi1)
d_psi12 = sym.simplify(d_psi1*d_psi2)

d_psi22 = sym.simplify(d_psi2*d_psi2)

b0 = sym.simplify(f*psi0)
b1 = sym.simplify(f*psi2)
b2 = sym.simplify(f*psi2)

print("(0,0) = ", sym.simplify(sym.integrate(d_psi00, (x, (1-(2*h)),(1-h)))))
print("(0,1) = ", sym.simplify(sym.integrate(d_psi01, (x, (1-(2*h)),(1-h)))))
print("(0,2) = ", sym.simplify(sym.integrate(d_psi02, (x, (1-(2*h)),1))))

print("(1,1) = ", sym.simplify(sym.integrate(d_psi11, (x, (1-(2*h)),(1-h)))))
print("(1,2) = ", sym.simplify(sym.integrate(d_psi12, (x, (1-(h)),1))))

print("(2,2) = ", sym.simplify(sym.integrate(d_psi22, (x, (1-h),1))))

print("b0 = ", sym.simplify(sym.integrate(b0, (x, 0, 1))), "- C(1-h)/2h*h - D*", sym.simplify(sym.integrate(d_psi0, (x,0,1))))
print("b1 = ", sym.simplify(sym.integrate(b1, (x, 0, 1))), "- C(2h-1)/h*h - D*", sym.simplify(sym.integrate(d_psi1, (x,0,1))))
print("b2 = ", sym.simplify(sym.integrate(b2, (x, 0, 1))), "- C(1-3h+2*h*h)/2h*h - D*", sym.simplify(sym.integrate(d_psi2, (x,0,1))))
