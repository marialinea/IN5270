import numpy as np
import matplotlib.pyplot as plt
import sympy as sym




class P2_FEM():
    def __init__(self, Ne, f, C, D, symbolic):
        self.Ne = Ne                    # Number of P2 elements
        self.f = lambda x: f(x)         # Right hand side of original PDE
        self.C = C                      # u'(0) = C
        self.D = D                      # u(1) = D
        self.step = 1/(2*Ne)            # step lenght between nodes
        self.Nn = 2*Ne+1                # number of nodes
        self.symbolic = symbolic

        self.nodes = np.linspace(0,1,self.Nn)
        self.element_matrix = np.zeros([2*self.Ne, 2*self.Ne])
        self.element_vector = np.array(2*Ne)

        self.psi_sym = sym.zeros(3*self.Ne,1)  # holds symbolic basis functions
        self.psi_num = sym.zeros(3*self.Ne,1)  # holds numeric basis functions
        self.d_psi_sym = sym.zeros(3*self.Ne,1)  # holds symbolic derivatives of basis functions
        self.d_psi_num = sym.zeros(3*self.Ne,1)  # holds numeric derivatives of basis functions

        self.x = sym.Symbol("x")
        self.h = sym.Symbol("h")
        self.c = sym.Symbol("C")
        self.d = sym.Symbol("D")



    def BasisFunctions(self):
        temp = -1
        for i in range(1, 2*self.Ne-2,2):
            self.psi_sym[temp + i] = (self.x-i*self.h)*(self.x-i*self.h-self.h)/(2*self.h*self.h)
            self.psi_sym[temp + i+1] = -(self.x-i*self.h+self.h)*(self.x-i*self.h-self.h)/(self.h*self.h)
            self.psi_sym[temp + i+2] = (self.x-i*self.h)*(self.x-i*self.h+self.h)/(2*self.h*self.h)

            self.d_psi_sym[temp + i] = sym.diff(self.psi_sym[i-1], self.x)
            self.d_psi_sym[temp + i+1] = sym.diff(self.psi_sym[i], self.x)
            self.d_psi_sym[temp + i+2] = sym.diff(self.psi_sym[i+1], self.x)

            temp += 1

        self.psi_sym[-3] = (self.x-1+self.h)*(self.x-1)/(2*self.h*self.h)
        self.psi_sym[-2] = -(self.x-1+2*self.h)*(self.x-1)/(self.h*self.h)
        self.psi_sym[-1] = (self.x-1+self.h)*(self.x-1+2*self.h)/(2*self.h*self.h)


        self.d_psi_sym[-3] = sym.diff(self.psi_sym[2*self.Ne-2], self.x)
        self.d_psi_sym[-2] = sym.diff(self.psi_sym[2*self.Ne-1], self.x)
        self.d_psi_sym[-1] = sym.diff(self.psi_sym[2*self.Ne], self.x)
        print(self.psi_sym)

        if not self.symbolic:
            self.psi_num = [sym.lambdify([self.x,self.h], self.psi_sym[i], modules="numpy") for i in range(self.Nn)]

            self.d_psi_num = [sym.lambdify([self.x,self.h], self.d_psi_sym[i], modules="numpy") for i in range(self.Nn)]
