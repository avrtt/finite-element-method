## General Introduction
<a href="https://www.youtube.com/watch?v=nxXcuDAv7Ss&t=107s"><img src="https://github.com/lenferdetroud/misc/blob/master/university/vordt_of_the_boreal_valley.png" alt="Vordt of the Boreal Valley" width="450"/></a>  

## Contents
1. [Problem Formulation](https://github.com/lenferdetroud/finite-element-method#problem-formulation)
2. [Introduction to Finite Element Method](https://github.com/lenferdetroud/finite-element-method#introduction-to-finite-element-method)
3. [Plan](https://github.com/lenferdetroud/finite-element-method#plan)
    * [Understanding the Problem](https://github.com/lenferdetroud/finite-element-method#understanding-the-problem)
    * [Galerkin Method](https://github.com/lenferdetroud/finite-element-method#galerkin-method)
    * [Basis Functions](https://github.com/lenferdetroud/finite-element-method#basis-functions)
    * [Local Matrices](https://github.com/lenferdetroud/finite-element-method#local-matrices)
    * [Stiffness and Mass Matrices](https://github.com/lenferdetroud/finite-element-method#stiffness-and-mass-matrices)
    * [Visualization](https://github.com/lenferdetroud/finite-element-method#visualization)
    * [Assembling the Global Matrix](https://github.com/lenferdetroud/finite-element-method#assembling-the-global-matrix)
    * [Boundary Conditions](https://github.com/lenferdetroud/finite-element-method#boundary-conditions)
4. [Implementation](https://github.com/lenferdetroud/finite-element-method#implementation)
5. [Testing](https://github.com/lenferdetroud/finite-element-method#testing)
6. [Research](https://github.com/lenferdetroud/finite-element-method#research)
7. [Results](https://github.com/lenferdetroud/finite-element-method#results)

## Problem Formulation
Implement finite element method for two-dimensional boundary value problem for elliptic equation in Cartesian coordinate system with bilinear basis functions on rectangles and boundary conditions of all types. Decompose the diffusion coefficient by biquadratic basis functions and generate the matrix of the system of linear algebraic equations in sparse string format. To solve the system of linear algebraic equations use the method of conjugate gradients or local-optimal scheme with incomplete factorization.
## Introduction to Finite Element Method
...
## Plan
#### Understanding the Problem
Equation of the elliptic boundary value problem:
[1.svg](https://github.com/lenferdetroud/misc/blob/master/finite-element-method/1.svg)
#### Galerkin Method
...
#### Basis Functions 
...
#### Local Matrices
...
#### Stiffness and Mass Matrices
...
#### Visualization
...
#### Assembling the Global Matrix
...
#### Boundary Conditions 
...
## Implementation
- [io.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/io.cpp)  
...
- [calculate\_local.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/calculate_local.cpp)  
...
- [generate\_global\_portrait.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/generate_global_portrait.cpp)  
...
- [assembling.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/assembling.cpp)  
...
- [boundary\_conditions.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/boundary_conditions.cpp)  
...
- [system\_solve.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/system_solve.cpp)  
...
- [main.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/main.cpp)  
...
- [data](https://github.com/lenferdetroud/finite-element-method/blob/main/)  
...
## Testing
...
## Research
...
## Results
...
