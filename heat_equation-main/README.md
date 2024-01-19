
# Steady state heat equation in a rectangle



## Intruduction

Hello, I am John Doe, and I am a student on the University of Science. This is my final project in the programming class.

It is a program which iteratively solves the steady state heat equation on a rectangular plate. The rectangle is dicretized into a grid of $nx \times ny$ points, in which the temperature of the plate is monitored.

On the east, west and south boundary, the rectangular plate is being heated to a temperature of 100 °C, while on the north boundary, the plate is being cooled to stay at 0 °C:
```
                   T = 0
             +------------------+
             |                  |
    T = 100  |                  | T = 100
             |                  |
             +------------------+
                   T = 100
```

The steady state solution to the discrete heat equation satisfies the condition, that in each interior point of the grid, the temperature is equal to the average temperature of the 4 neighboring grid points:
```
T[central] = (1/4) * ( T[north] + T[south] + T[east] + T[west] )
```

Given an approximate solution of the steady state heat equation, a better solution can be obtained by replacing each interior point by the average of its 4 neighbors - that is, by using the above condition as an assignment. If we iteratively repeat this process enough times, we can reach the desired convergence criteria, which is satisfied if the maximum difference between temperatures in successive iterations reaches the desired tolerance.



## Description of the project

The main program in this project is `heat_equation`, which performs the iteration until convergence is reached. It creates the grid of points, initializes it with the boundary conditions and an initial estimate of the solution, solves the heat equation using the iterative method, and writes the solution to a given binary file. The other program, `heat_to_bmp`, which I copied from somewhere, converts the solution to a bitmap image.

To make the programs work, I need to execute this command first
```
module load GCC/13.2.0
```

Create a directory for the output files
```
mkdir io
```

To compile the programs, I use
```
g++ -g -O2 src/heat_equation.cpp -o heat_equation
g++ -g -O2 src/heat_to_bmp.cpp -o heat_to_bmp
```

To solve the discrete steady-state heat equation on a grid of 1200-by-1000 ($nx$-by-$ny$) points, use e.g.
```
./heat_equation 1200 1000 io/heat.bin
```
It takes just over a minute to run this program.

To then convert the solution to a bitmap image, use
```
./heat_to_bmp io/heat.bin io/heat.bmp
```



## Task

Unfortunately, my programming teacher was not happy with this project. He said that the conjugate gradient solver is very slow. But I have no idea how to make it faster.

__Will you help me make the program as fast as possible?__

You can modify my code as much as you like. But the main functionality of the program has to stay the same - iteratively solve the discrete steady state heat equation on a recrangular grid.
