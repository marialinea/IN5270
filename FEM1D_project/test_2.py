import matplotlib.pyplot as plt
import numpy as np

def analytical(x, C, D):
    return 0.5*x*x - 1/3 * x**3 + C*x + D - C - 1/6

N = 2
C = 0
D = 0
h = 1/(2*N)

A = np.array([[13/12, -7/6, 1/6, 0,0,0], [-7/6, 4/3, -7/6,0,0,0], [1/6, -7/6, 13/12,0,0,0], [0,0,0,13/12, -7/6, 1/6], [0,0,0,-7/6, 4/3, -7/6], [0,0,0,1/6, -7/6, 13/12]]).reshape(6,6) * 4

b = np.array([3/(4*h) - (6*h+1)/(6*h*h) + 1/(4*h*h), 1/(4*h) - (2*h+1)/(6*h*h) + 1/(4*h*h), 1/(4*h) - (2*h+1)/(6*h*h) + 1/(4*h*h), (h-1)/(12*h*h), (3*h-1)/(12*h*h), (3*h-1)/(12*h*h)]).reshape(6,1)

print(b)

X = np.linalg.solve(A,b)
print(X)



def u(x):
    psi0 = (x-h)*(x-2*h)/(2*h*h)
    psi1 = -x*(x-2*h)/(h*h)
    psi2 = x*(x-h)/(2*h*h)

    psiNm2 = (x-1+h)*(x-1)/(2*h*h)
    psiNm1 = -(x-1+2*h)*(x-1)/(h*h)
    psiN = (x-1+2*h)*(x-1+h)/(2*h*h)
    if x <= 0.5:
        return X[0]*psi0 + X[1]*psi1 + X[2]*psi2
    else:
        return X[3]*psiNm2 + X[4]*psiNm1 + X[5]*psiN


u = np.vectorize(u)

x = np.linspace(0,1,1000)
plt.plot(x, analytical(x, C, D), label="exact")
plt.plot(x, u(x), label="numerical")
plt.legend()
plt.show()
