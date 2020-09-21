import numpy as np
import matplotlib.pyplot as plt

class OneDim_WaveSolver():

    def __init__(self, Nx, L, c, T):
        self.Nx = Nx
        self.L = L
        self.c = c
        self.T = T

        self.u_new = np.zeros(Nx)
        self.u = np.zeros(Nx)
        self.u_old = np.zeros(Nx)

        self.dx = L/(Nx)
        self.C = 0.5
        self.dt = self.C*self.dx/c
        self.Nt = T/self.dt
        self.x = np.linspace(0,L,Nx)

        self.l2_norm = 0
        self.l_infty_norm = 0

    def Conditions(self, I):
        self.u_old = I(self.x)

        #Boundary conditions
        self.u[0] = 0
        self.u[-1] = 0

        self.u_old[0] = 0.
        self.u_old[-1] = 0.

        self.u_new[0] = 0

    def FirstTimeStep(self):
        for i in range(1,self.Nx-1):
            self.u[i] = self.u_old[i] - 0.5*self.C**2*(self.u_old[i+1] - 2*self.u_old[i] + self.u_old[i-1])

    def AdvanceTime(self, exact):
        t = 2*self.dt
        while t <= self.T:
            self.AdvanceSpace()
            self.VectorSwap()
            self.L2norm(exact, t)
            t += self.dt
        self.l2_norm = np.sqrt(self.dt*self.dx*self.l2_norm)

    def VectorSwap(self):
        self.u_old[:] = self.u[:]
        self.u[:] = self.u_new[:]

        """
        Raskere versjon av swapping:
        tmp = u_old
        u_old = u
        u = u_new
        u_new = tmp     # skriver over alle disse verdiene fordi de er ubrukelige
        # u_new skal uansett oppdateres, men det går minst tid å bare sette det til u_old.

        """

    def AdvanceSpace(self):
        u_new, u, u_old, C = self.u_new, self.u, self.u_old, self.C
        for i in range(1,self.Nx-1):
            u_new[i] = -u_old[i] + 2*u[i] + C**2 * (u[i+1] - 2*u[i] + u[i-1])


        """
        Vektorisert versjon:
        u_new[1:-1] = -u_old[1:-1] + 2*u[1:-1] + C**2 * (u[2:] - 2*u[1:-1] + u[:-2])
        """


    def L2norm(self, exact, t):
        difference = exact(self.x, t) - self.u[:]
        self.l2_norm += np.sum(difference**2)




    def PlotComparing(self, exact):
        plt.plot(self.x, exact(self.x, self.T), label="Exact")
        plt.plot(self.x, self.u, label="Numerical")
        plt.legend()
        plt.xlabel("x")
        plt.ylabel("u(x,T)")
        plt.show()
