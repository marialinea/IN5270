import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy.printing.latex import LatexPrinter, print_latex



class P2_FEM():
    def __init__(self, Ne, f, C, D, analytical):
        """
        Input parameters:
        -----------------

        Ne : int
            Number of P2 elements
        f : string
            The right-hand side in eq. -u''(x) = f(x), given as a string
        C : int
            Initial condition u'(0) = C
        D : int
            Boundary condition u(1) = D
        analytical : string
            Analytical solution given as a string, as a function of x, C and D

        """
        self.Ne = Ne                    # Number of P2 elements
        self.f = sym.sympify(f)         # Right hand side of original PDE
        self.C = C                      # u'(0) = C
        self.D = D                      # u(1) = D
        self.step = 1/(2*Ne)            # step lenght between nodes

        self.element_matrix = sym.zeros(2*self.Ne, 2*self.Ne)
        self.element_vector = sym.zeros(2*Ne,1)

        self.psi_sym = sym.zeros(3*self.Ne,1)    # holds symbolic basis functions
        self.psi_num = sym.zeros(3*self.Ne,1)    # holds numeric basis functions
        self.d_psi_sym = sym.zeros(3*self.Ne,1)  # holds symbolic derivatives of basis functions
        self.d_psi_num = sym.zeros(3*self.Ne,1)  # holds numeric derivatives of basis functions

        self.x = sym.Symbol("x")
        self.h = sym.Symbol("h")
        self.c = sym.Symbol("C")
        self.d = sym.Symbol("D")

        self.analytical = sym.lambdify([self. x, self.c, self.d], analytical, modules="numpy")




    def DOFmap(self):
        self.dof_map = []

        node = 0
        for e in range(self.Ne):
            self.dof_map.append([])
            for i in range(node, 3+node):
                self.dof_map[e].append(i)
            node += 2

    def BasisFunctions(self):

        temp = -1
        for i in range(1, 2*self.Ne-2,2):
            self.psi_sym[temp + i] = (self.x-i*self.h)*(self.x-i*self.h-self.h)/(2*self.h*self.h)
            self.psi_sym[temp + i+1] = -(self.x-i*self.h+self.h)*(self.x-i*self.h-self.h)/(self.h*self.h)
            self.psi_sym[temp + i+2] = (self.x-i*self.h)*(self.x-i*self.h+self.h)/(2*self.h*self.h)

            self.d_psi_sym[temp + i] = sym.diff(self.psi_sym[temp + i], self.x)
            self.d_psi_sym[temp + i+1] = sym.diff(self.psi_sym[temp + i+1], self.x)
            self.d_psi_sym[temp + i+2] = sym.diff(self.psi_sym[temp + i+2], self.x)

            temp += 1

        self.psi_sym[-3] = (self.x-1+self.h)*(self.x-1)/(2*self.h*self.h)
        self.psi_sym[-2] = -(self.x-1+2*self.h)*(self.x-1)/(self.h*self.h)
        self.psi_sym[-1] = (self.x-1+self.h)*(self.x-1+2*self.h)/(2*self.h*self.h)


        self.d_psi_sym[-3] = sym.diff(self.psi_sym[-3], self.x)
        self.d_psi_sym[-2] = sym.diff(self.psi_sym[-2], self.x)
        self.d_psi_sym[-1] = sym.diff(self.psi_sym[-1], self.x)


        self.psi_num = [sym.lambdify([self.x,self.h], self.psi_sym[i], modules="numpy") for i in range(3*self.Ne)]
        self.d_psi_num = [sym.lambdify([self.x,self.h], self.d_psi_sym[i], modules="numpy") for i in range(3*self.Ne)]



    def MatrixAssembly(self):

        self.DOFmap()
        d_psi = self.d_psi_sym

        counter = 0
        for e in range(self.Ne):
            for i in self.dof_map[e]:
                if i == self.Ne*2:
                    break
                for j in range(i, self.dof_map[e][-1]+1):
                   if j == 2*self.Ne:
                       continue
                   self.element_matrix[i,j] += sym.integrate(sym.simplify(d_psi[i+counter]*d_psi[j+counter]), (self.x, self.h*self.dof_map[e][0],self.h*self.dof_map[e][-1]))
                   self.element_matrix[j,i] = self.element_matrix[i,j]
            counter += 1

        print_latex(self.element_matrix)
        self.element_matrix = np.array(self.element_matrix.subs(self.h,self.step)).astype(np.float64)

    def ElementVector(self):
        psi = self.psi_sym
        # Filling the first element
        self.element_vector[0] = sym.simplify(sym.integrate(self.f*psi[0], (self.x, 0, 2*self.h)) - self.c)

        # Filling the entries from the interior elements
        counter = 0
        for e in range(self.Ne-1):
            for i in self.dof_map[e]:
                if i == 0: continue
                if i == 2*self.Ne: break
                self.element_vector[i] += sym.simplify(sym.integrate(self.f*psi[i+counter], (self.x, self.h*self.dof_map[e][0],self.h*self.dof_map[e][-1])))
            counter += 1

        # Filling the entries from the last element

        self.element_vector[-2] += sym.simplify(sym.integrate(self.f*psi[-3], (self.x, self.h*self.dof_map[-1][0],self.h*self.dof_map[-1][-1])) - self.d*sym.integrate(self.d_psi_sym[-1]*self.d_psi_sym[-3],  (self.x, self.h*self.dof_map[-1][0],self.h*self.dof_map[-1][-1])))


        self.element_vector[-1] += sym.simplify(sym.integrate(self.f*psi[-2], (self.x, self.h*self.dof_map[-1][0],self.h*self.dof_map[-1][-1])) - self.d*sym.integrate(self.d_psi_sym[-1]*self.d_psi_sym[-2],  (self.x, self.h*self.dof_map[-1][0],self.h*self.dof_map[-1][-1])))

        self.element_vector = np.array(self.element_vector.subs({self.h:self.step, self.c:self.C, self.d:self.D})).astype(np.float64)



    def SolveLinearSystem(self):

        self.BasisFunctions()
        self.MatrixAssembly()
        self.ElementVector()

        self.weights = np.linalg.solve(self.element_matrix, self.element_vector)


        self.N = 1000
        self.solution = np.zeros(self.N*self.Ne)

        counter = 0
        for index, e in enumerate(range(self.Ne)):
            x = np.linspace(self.step*self.dof_map[e][0],self.step*self.dof_map[e][-1],self.N)

            term1 = self.weights[self.dof_map[e][0]]*self.psi_num[self.dof_map[e][0]+counter](x,self.step)
            term2 = self.weights[self.dof_map[e][1]]*self.psi_num[self.dof_map[e][1]+counter](x,self.step)

            if e == self.Ne-1:
                term3 = self.psi_num[self.dof_map[e][2]+counter](x,self.step)*self.D
            else:
                term3 = self.weights[self.dof_map[e][2]]*self.psi_num[self.dof_map[e][2]+counter](x,self.step)

            self.solution[index*self.N:(index+1)*self.N] = term1 + term2 + term3
            counter += 1

    def PlotSolution(self):

        x = np.linspace(0,1,self.N*self.Ne)
        plt.plot(x, self.solution, label="Approximation")
        plt.plot(x, self.analytical(x, self.C, self.D), label="Analytical")
        plt.xlabel("x"); plt.ylabel("u(x)")
        plt.legend()
        plt.show()


    def L2Norm(self):

        x = np.linspace(0,1,self.N*self.Ne)
        self.L2 = np.linalg.norm((self.analytical(x, self. C, self.D)-self.solution))
