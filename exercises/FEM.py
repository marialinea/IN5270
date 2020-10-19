import matplotlib.pyplot as plt
import numpy as np

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
