import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import fractions

x = sym.Symbol("x")
h = sym.Symbol("h")
f = sym.sympify("1 + 2*x - x*x")

def Lagrange_polynomial(x, i, points):
    p=1
    for k in range(len(points)):
        if k != i:
            p *= (x - points[k])/(points[i] - points[k])
    return p

def DOFmap(Num_elem, deg):
    Ne = Num_elem
    dof_map = []

    node = 0
    if Ne == 1:
        for i in range(node, deg+1+node):
            dof_map.append(i)
    else:
        for e in range(Ne):
            dof_map.append([])
            for i in range(node, deg+1+node):
                dof_map[e].append(i)
            node += deg
    return dof_map


element_matrix = np.zeros([3,3])
element_vector = np.zeros(3)

# One P2 element
"""

dof_map = DOFmap(1, 2)

points = [0,0.5,1]

psi = [Lagrange_polynomial(x, r, points) for r in range(3)]

counter = 0
for i in dof_map:
    for j in range(i, dof_map[-1]+1):
       element_matrix[i,j] += sym.integrate(sym.simplify(psi[i+counter]*psi[j+counter]), (x,0,1))
       element_matrix[j,i] = element_matrix[i,j]
counter += 1


# Filling the entries from the interior elements
counter = 0
for i in dof_map:
    element_vector[i] += sym.simplify(sym.integrate(f*psi[i+counter], (x,0,1)))
counter += 1


u = np.linalg.solve(element_matrix, element_vector)

psi_num = [sym.lambdify([x,h], psi[i], modules="numpy") for i in range(3)]
h = 1/2


def f(x):
    return 1 + 2*x - x*x
x = np.linspace(0,1,1000)

plt.plot(x, f(x), label="analytical")
plt.plot(x, u[0]*psi_num[0](x,h) + u[1]*psi_num[1](x,h) + u[2]*psi_num[2](x,h), label="FEM")
plt.legend()
plt.show()
"""

# Two P1 elements

psi = sym.zeros(4,1)

x1 = [0, 0.5]
x2 = [0.5, 1]

psi[0] = Lagrange_polynomial(x, 0, x1)
psi[1] = Lagrange_polynomial(x, 1, x1)
psi[2] = Lagrange_polynomial(x, 0, x2)
psi[3] = Lagrange_polynomial(x, 1, x2)


dof_map = DOFmap(2, 1)


counter = 0
Ne = 2
step = 1/Ne
for e in range(Ne):
    for i in dof_map[e]:
        for j in range(i, dof_map[e][-1]+1):
           element_matrix[i,j] += sym.integrate(sym.simplify(psi[i+counter]*psi[j+counter]), (x, step*dof_map[e][0],step*dof_map[e][-1]))
           element_matrix[j,i] = element_matrix[i,j]
    counter += 1


element_matrix = np.array(element_matrix).astype(np.float64)



# Filling the entries from the interior elements
counter = 0
for e in range(Ne):
    for i in dof_map[e]:
        element_vector[i] += sym.simplify(sym.integrate(f*psi[i+counter], (x, step*dof_map[e][0],step*dof_map[e][-1])))
    counter += 1

element_vector = np.array(element_vector).astype(np.float64)


psi_num = [sym.lambdify([x,h], psi[i], modules="numpy") for i in range(4)]

def f(x):
    return 1 + 2*x - x*x

def u(x):
    u = np.linalg.solve(element_matrix, element_vector)

    psi1 = psi_num[0](x,step)
    psi2 = psi_num[1](x,step)
    psi3 = psi_num[2](x,step)
    psi4 = psi_num[3](x,step)
    if x <= 0.5:
        return u[0]*psi1 + u[1]*psi2
    else:
         return u[1]*psi3 + u[2]*psi4


u = np.vectorize(u)

x = np.linspace(0,1,1000)
plt.plot(x, f(x), label="analytical")
plt.plot(x, u(x), label ="FEM")
plt.show()
