import matplotlib.pyplot as plt
import numpy as np

def analytical(x, C, D):
    return 0.5*x*x - 1/3 * x**3 + C*x +D -C -1/6

N = 1
C = 0
D = 0
h = 1/(2*N)

A = np.array([[13/12, -7/6, 1/6], [-7/6, 4/3, -7/6], [1/6, -7/6, 13/12]]).reshape(3,3) * 2

b = np.array([3/(4*h) - (6*h+1)/(6*h*h) + 1/(4*h*h), 1/(4*h) - (2*h+1)/(6*h*h) + 1/(4*h*h), 1/(4*h) - (2*h+1)/(6*h*h) + 1/(4*h*h)]).reshape(3,1)

X = np.linalg.solve(A,b)
print(X)



def u(x):
    """
    psi0 = (x-1+h)*(x-1)/(2*h*h)
    psi1 = -(x-1+2*h)*(x-1)/(h*h)
    psi2 = (x-1+2*h)*(x-1+2*h)/(2*h*h)
    """
    psi0 = (x-h)*(x-2*h)/(2*h*h)
    psi1 = -x*(x-2*h)/(h*h)
    psi2 = x*(x-h)/(2*h*h)
    return X[0]*psi0 + X[1]*psi1 + X[2]*psi2


u = np.vectorize(u)

x = np.linspace(0,1,1000)
plt.plot(x, analytical(x, C, D), label="exact")
plt.plot(x, u(x), label="numerical")
plt.legend()
plt.show()
