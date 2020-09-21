import numpy as np
import mayavi.mlab as mlab



class WaveSolver2D():

    def __init__(self, b, Nx, Ny, Lx, Ly, T):
        self.b = b
        self.Nx = Nx
        self.Ny = Ny
        self.Lx = Lx
        self.Ly = Ly
        self.T = T

        self.x = np.linspace(0,self.Lx,self.Nx)
        self.y = np.linspace(0,self.Ly,self.Ny)
        self.dx = self.x[1]-self.x[0]
        self.dy = self.y[1]-self.y[0]


        self.u = np.zeros([self.Nx+2, self.Ny+2])
        self.u_old = np.zeros([self.Nx+2, self.Ny+2])
        self.u_new = np.zeros([self.Nx+2, self.Ny+2])

    def Initialize(self, I, V, q, f):
        self.I = lambda x,y: I(x,y)
        self.V = lambda x,y: V(x,y)
        self.q_func = lambda x,y: q(x,y)
        self.f = lambda x,y,t: f(x,y,t)

        self.dt = 0.9* (1/np.sqrt(np.max(self.q_func(self.x,self.y)) *(1/self.dx**2 + 1/self.dy**2)))
        self.Nt = self.T/self.dt
        self.dtdt = self.dt**2
        self.dxdx = self.dx**2
        self.dydy = self.dy**2
        self.A = 1/(2+self.b*self.dt)

        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.X = self.X.T; self.Y = self.Y.T

        self.SetInitialCondition()
        self.MakeQMatrix()


    def SetInitialCondition(self):
        """
        Setting initial condition for the spatial coordinates.
        """

        # Filling the interior, including the boundary, points
        for i in range(1,self.Nx+1):
            for j in range(1,self.Ny+1):
                self.u_old[i,j] = self.I(self.x[i-1],self.y[j-1])

        self.UpdateGhostCells(self.u_old)


    def MakeQMatrix(self):
        """
        Stores the function values of q in a matrix. Reduces the number of times
        we call the function during our calculation.
        """

        # Filling the interior, including the boundary, points
        self.q = np.zeros([self.Nx+2, self.Ny+2])
        for i in range(1,self.Nx+1):
            for j in range(1,self.Ny+1):
                self.q[i,j] = self.q_func(self.x[i-1],self.y[j-1])

        # Setting up the ghost values
        for i in range(1,self.Nx+1):
            self.q[i,0] = 2*self.q[i,1]-self.q[i,2]
            self.q[i,self.Ny+1] = 2*self.q[i,self.Ny]-self.q[i,self.Ny-1]

        for j in range(1,self.Ny+1):
            self.q[0,j] = 2*self.q[1,j]-self.q[2,j]
            self.q[self.Nx+1,j] = 2*self.q[self.Nx,j]-self.q[self.Nx-1,j]

    def UpdateGhostCells(self, u):
        """
        Updates the ghost values.
        """

        for i in range(1,self.Nx+1):
            u[i,0] = u[i,2]
            u[i,self.Ny+1] = u[i,self.Ny-1]

        for j in range(1,self.Ny+1):
            u[0,j] = u[2,j]
            u[self.Nx+1,j] = u[self.Nx-1,j]

    def UpdateGhostCells_vec(self, u):

        u[1:-1,0] = u[1:-1,2]
        u[1:-1,self.Ny+1] = u[1:-1,self.Ny-1]

        u[0,1:-1] = u[2,1:-1]
        u[self.Nx+1,1:-1] = u[self.Nx-1,1:-1]

    def FirstTimeStep(self):
        u_old = self.u_old
        q = self.q

        for i in range(1,self.Nx+1):
            for j in range(1, self.Ny+1):
                self.u[i,j] = 0.25*(2*(2-self.b*self.dt)*self.dt*self.V(self.x[i-1],self.y[j-1]) + (self.dt/self.dx)**2 * ((q[i+1,j]+q[i,j])*(u_old[i+1,j]-u_old[i,j])-(q[i,j]+q[i-1,j])*(u_old[i,j]-u_old[i-1,j])) + (self.dt/self.dy)**2 * ((q[i,j+1]+q[i,j])*(u_old[i,j+1]-u_old[i,j])-(q[i,j]+q[i,j-1])*(u_old[i,j]-u_old[i,j-1])) + 2 * self.dt**2 * self.f(self.x[i-1],self.y[j-1],0.) + 4*u_old[i,j] )


        self.UpdateGhostCells_vec(self.u)



    def AdvanceTime(self):
        u_old = self.u_old
        q = self.q

        counter = 0
        self.t = self.dt
        while self.t <= self.T:
            self.AdvanceSpace_vec(self.t)
            if counter%1000 == 0:
                filename = "{}.png".format(counter)
                self.Animate(filename)
            self.VectorSwap()
            self.t += self.dt
            counter += 1


    def AdvanceSpace(self, t):
        u, u_old, q = self.u, self.u_old, self.q

        for i in range(1, self.Nx+1):
            for j in range(1, self.Ny+1):
                self.u_new[i,j] = self.A*(self.b*self.dt - 2)*u_old[i,j] + 4*self.A*u[i,j] + (self.dt/self.dx)**2 * self.A * ((q[i+1,j]+q[i,j])*(u[i+1,j]-u[i,j])-(q[i,j]+q[i-1,j])*(u[i,j]-u[i-1,j])) + (self.dt/self.dy)**2 * self.A * ((q[i,j+1]+q[i,j])*(u[i,j+1]-u[i,j])-(q[i,j]+q[i,j-1])*(u[i,j]-u[i,j-1])) + 2*self.A*self.dt**2 * self.f(self.x[i-1], self.y[j-1], t)

        self.UpdateGhostCells_vec(self.u_new)

    def AdvanceSpace_vec(self, t):
        u, u_old, q = self.u, self.u_old, self.q

        self.u_new[1:-1,1:-1] = self.A*(self.b*self.dt - 2)*u_old[1:-1,1:-1] + 4*self.A*u[1:-1,1:-1] + (self.dt/self.dx)**2 * self.A * ((q[2:,1:-1]+q[1:-1,1:-1])*(u[2:,1:-1]-u[1:-1,1:-1])-(q[1:-1,1:-1]+q[:-2,1:-1])*(u[1:-1,1:-1]-u[:-2,1:-1])) + (self.dt/self.dy)**2 * self.A * ((q[1:-1,2:]+q[1:-1,1:-1])*(u[1:-1,2:]-u[1:-1,1:-1])-(q[1:-1,1:-1]+q[1:-1,:-2])*(u[1:-1,1:-1]-u[1:-1,:-2])) + 2*self.A*self.dt**2 * self.f(self.X, self.Y, t)

        self.UpdateGhostCells_vec(self.u_new)

    def VectorSwap(self):
        tmp = self.u_old
        self.u_old = self.u
        self.u = self.u_new
        self.u_new = tmp

    def ComputeError(self, exact):
        self.exact = lambda x,y,t: exact(x,y,t)

        exact_sol = np.zeros([self.Nx, self.Ny])
        for i in range(self.Nx):
            for j in range(self.Ny):
                exact_sol[i,j] = self.exact(self.x[i],self.y[j],self.t)

        l_infty_norm = exact_sol - self.u[1:-1,1:-1]
        return np.max(np.abs(l_infty_norm))


    def Animate(self, filename):
        mlab.clf()
        extent1 = (0, 20, 0, 20,0, 10)
        s = mlab.surf(self.x , self.y, self.u[1:-1,1:-1], colormap="Blues", warp_scale=5,extent=extent1)
        mlab.axes(s, color=(.7, .7, .7), extent=extent1, ranges=(0, 20, 0, 20, 0, 10), xlabel="", ylabel="", zlabel="", x_axis_visibility=False, z_axis_visibility=False)
        mlab.outline(s, color=(0.7, .7, .7), extent=extent1)
        mlab.text(6, -2.5, "", z=-4, width=0.14)
        mlab.colorbar(object=None, title=None, orientation="horizontal", nb_labels=None, nb_colors=None, label_fmt=None)
        mlab.title("Gaussian t=%g" % self.t)
        mlab.view(142, -72, 50)
        f = mlab.gcf()
        camera = f.scene.camera
        camera.yaw(0)
        mlab.savefig(filename)
