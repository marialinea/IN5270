## Short overview of the files used in this project

All of the codes are found [here](https://github.com/marialinea/IN5270/tree/master/FEM1D_project).

* FEM_class.py

 This file contains a solver specialized for P2 elements for a PDE with an analytical solution on the domain x in [0,1].

 It contains the follwing methods:
  1. DOFmap

    Constructs a degrees of freedom map. Holds the global node number for each element.
  2. BasisFunctions

    Constructs second order Lagrange polynomials for each node.

  3. MatrixAssembly

    Assembles the coefficient matrix, A.

  4. ElementVector

    Piecing together the element vector, b.

  5. SolveLinearSystem

    Solves the linear system Ac = b, where c is the coefficient vector that contains the coefficients in the linear combination for u.

  6. PlotSolution

      Plots the numerical approximation vs. the analytical solution.

  7. L2Norm

    Computes the L2 norm.

* main.py

 This program calls the class, and solves the different tasks in the project.
