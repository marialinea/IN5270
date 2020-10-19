## Short overview of the files used in this project

All of the codes are found [here](https://github.com/marialinea/IN5270/tree/master/Project%201%20-%20Wave%20Equation/codes).

* WaveSolver2D.py

 This file contains the general 2D wave equation solver. It advances the wave equation in both time and space, and the class method that evolves the equation in space has a vectorized version.

 The final solution is for the final time step, i. e. the solver do not store the solution u for all time steps.

 The solver also generates the png files used to make the gifs in part 4 of the project.

* animate.py

 Concatenates all the png files into a gif and deletes all the png files.

* main.py

 This is the main program that calls the class WaveSolver2D, and solves the different tasks in the project. To see each code segment that are used in the different parts of the project see the final report [here](https://github.com/marialinea/IN5270/blob/master/Project%201%20-%20Wave%20Equation/report/Project%201%20-%20Finite%20difference%20simulation%20of%202D%20waves.pdf).
