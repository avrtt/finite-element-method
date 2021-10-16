## General Introduction
<a href="https://www.youtube.com/watch?v=nxXcuDAv7Ss&t=107s"><img src="https://github.com/lenferdetroud/misc/blob/master/university/vordt_of_the_boreal_valley.png" alt="Vordt of the Boreal Valley" width="450"/></a>  
## Contents
...
## Problem Formulation
Implement finite element method for two-dimensional boundary value problem for elliptic equation in Cartesian coordinate system with bilinear basis functions on rectangles and boundary conditions of all types. Decompose the diffusion coefficient by biquadratic basis functions and generate the matrix of the system of linear algebraic equations in sparse string format. To solve the system of linear algebraic equations use the method of conjugate gradients or local-optimal scheme with incomplete factorization.
## Introduction to Finite Element Method
...
## Plan
#### Understanding the Problem
  
#### Galerkin Method
  
#### Basis Functions 
  
#### Local Matrices
  
#### Stiffness and Mass Matrices
  
#### Visualization
  
#### Assembling the Global Matrix
  
#### Boundary Conditions 
  
## Implementation
#### [io.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/io.cpp)
...
#### [calculate\_local.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/calculate_local.cpp)
...
#### [generate\_global_\portrait.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/generate_global_portrait.cpp)
...
#### [assembling.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/assembling.cpp)
...
#### [boundary\_conditions.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/boundary_conditions.cpp)
...
#### [system\_solve.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/system_solve.cpp)
...
#### [main.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/main.cpp)
...
#### [data](https://github.com/lenferdetroud/finite-element-method/blob/main/)
...
## Testing
...
## Research
...
## Results
...
## Literature
...
