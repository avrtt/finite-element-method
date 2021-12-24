<img src="https://github.com/lenferdetroud/misc/blob/master/fem/mesh.png" alt="mesh" align="right" width="370" height="305">

This repository contains an implementation of the two-dimensional [finite element method](https://en.wikipedia.org/wiki/Finite_element_method) in the Cartesian coordinate system.  
  
The FEM is a numerical method that solves [boundary value problems](https://en.wikipedia.org/wiki/Boundary_value_problem) (i.e., it solves [partial differential equations](https://en.wikipedia.org/wiki/Partial_differential_equation)) by discretizing the area. This is achieved by the [construction of a mesh](https://en.wikipedia.org/wiki/Mesh_generation) of the object: the numerical domain for the solution, which has a finite number of points. The FEM formulation of a boundary value problem finally results in a system of algebraic equations. The simple equations that model these finite elements are then assembled into a larger system of equations that models the entire problem. The method approximates the unknown function over the domain.


## Contents
1. [Outline](https://github.com/lenferdetroud/finite-element-method#outline)
    * [Understanding the Problem](https://github.com/lenferdetroud/finite-element-method#understanding-the-problem)
    * [Galerkin Method](https://github.com/lenferdetroud/finite-element-method#galerkin-method)
    * [Basis Functions](https://github.com/lenferdetroud/finite-element-method#basis-functions)
    * [Local Matrices](https://github.com/lenferdetroud/finite-element-method#local-matrices)
    * [Global Matrix](https://github.com/lenferdetroud/finite-element-method#global-matrix)
    * [Sparse String Format](https://github.com/lenferdetroud/finite-element-method#sparse-string-format)
    * [Boundary Conditions](https://github.com/lenferdetroud/finite-element-method#boundary-conditions)
    * [LU Decomposition and Gaussian Elimination](https://github.com/lenferdetroud/finite-element-method#lu-decomposition-and-gaussian-elimination)
    * [Local Optimal Scheme](https://github.com/lenferdetroud/finite-element-method#local-optimal-scheme)
2. [Implementation](https://github.com/lenferdetroud/finite-element-method#implementation)
3. [Testing and Research](https://github.com/lenferdetroud/finite-element-method#testing-and-research)
4. [Install](https://github.com/lenferdetroud/finite-element-method#install)


## Outline
### Understanding the Problem
We will consider a two-dimensional [elliptic equation](https://en.wikipedia.org/wiki/Elliptic_partial_differential_equation):   
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/9.png" width="35%">  
Using the definitions of [gradient](https://en.wikipedia.org/wiki/Gradient#Cartesian_coordinates) and [divergence](https://en.wikipedia.org/wiki/Divergence#Cartesian_coordinates) (in the case of two variables in the Cartesian coordinate system), the equation can be rewritten as follows:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/11.png" width="35%">  
This equation is set in some area **Ω** with boundary  <img src="https://github.com/lenferdetroud/misc/blob/master/fem/10.png" width="14%">  and has boundary conditions:   
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/12.png" width="60%">  
  
Our goal is to find **u = u(x, y)** (i.e., its numerical approximation).
The first step is to find an arbitrary trial function which satisfies the given boundary conditions. Substituting it into the equation, we will calculate the [residual](https://en.wikipedia.org/wiki/Residual_(numerical_analysis)) (the difference between the chosen and theoretical functions) to estimate the accuracy of the approximation. Next, we need to find such a function among some set of functions that approximates the **u** in the best way, i.e. solve the problem of minimization of a [functional](https://en.wikipedia.org/wiki/Functional_(mathematics)). 

### Galerkin Method
The best method for minimizing the error in our case is the [Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method), which can be described by this equation:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/0.png" width="20%">  
where **R0** is an error (residual) and **Ψ** is some arbitrary function **coinciding with the trial function**.  

Let's take the left side of our elliptic equation as **R0** and divide the integral by the sum of integrals:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/18.png" width="55%">  
The FEM proceeds from a different form of the DE, called the [weak formulation](https://en.wikipedia.org/wiki/Weak_formulation). In order to get it, we will use [integration by parts](https://en.wikipedia.org/wiki/Integration_by_parts), and also take the boundary conditions into account. Let's use integration by parts in the first integral:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/19.png" width="45%">  
The integral on the boundary **S** is divided into three boundaries (**S1**, **S2**, **S3**) on which we have boundary conditions:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/43.png" width="20%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/44.png" width="30%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/45.png" width="40%">  
And then:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/20.png" width="85%">   
This is where the division into finite elements begins. We can now represent the function **u** as a decomposition of the basis functions **ψ** with weights **q**:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/21.png" width="15%">   
Now the approximation of the function **u** will be performed by docked local basis functions on finite elements, and the solution of the problem is a vector of weights (**q1**, **q2**, ... **qn**), which can be obtained by solving a system of linear equations.  
Let's substitute the decomposition of the function **u** into the Galerkin equation and obtain the final expression for the linear system with boundary conditions:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/42.png" width="100%">  

### Basis Functions 
We will represent the finite elements as **rectangles**, and as functions **Ψ** we will take **bilinear functions**. The area will be divided into subareas corresponding to each finite element:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/22.png" width="10%">  
Bilinear functions are the product of the following one-dimensional functions on the corresponding axes:  <img src="https://github.com/lenferdetroud/misc/blob/master/fem/17.png" align="right" width="30%">   
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/15.png" width="63%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/16.png" width="60%">  
So, there are four local basis functions (on one finite element):   
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/14.png" width="50%">  
  
Each local basis function equals **1** at one node of its finite element **Ωk**.  
At other nodes and finite elements it equals **0**.  <img src="https://github.com/lenferdetroud/misc/blob/master/fem/23.png" align="right" width="30%">   
    
In addition to the bilinear functions, we similarly construct **biquadratic functions** on each finite element, which we will use to decompose the diffusion coefficient **λ**:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/24.png" width="60%">  
  
The corresponding one-dimensional functions:  
  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/71.png" width="90%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/75.png" width="90%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/77.png" width="75%">  


### Local Matrices
Now we need to proceed to the solution on each finite element from which the solution of the whole problem will be assembled. For this, let's present the equation obtained by Galerkin's method as a sum of integrals over regions **Ωk** without taking into account boundary conditions, and let's also disclose the gradients:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/39.png" width="65%">  
The first term gives the [stiffness matrix](https://en.wikipedia.org/wiki/Stiffness_matrix) **G**, the second term gives the [mass matrix](https://en.wikipedia.org/wiki/Mass_matrix) **M**. The sum of these matrices is called a local matrix **A=G+M** and has dimension **4x4** (based on the number of nodes). The integral on the right side is called a **local vector b**.   
  
It's time to build **G**, **M**, and **b**. Let's decompose the diffusion coefficient (i.e., construct an interpolating function) and derive the formula for calculating the elements of the matrix **G**:   
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/27.png" width="13%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/30.png" width="65%">  
In this expression **λ1**, **λ2**, ..., **λ9** are the values of **λ** at the corresponding nodes of the finite element.  
The matrix is too large to be shown here because each of its 16 elements has 9 terms.  
  
For the matrix **M**, we replace the parameter **γ** by its mean value on each finite element **Ωk**:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/37.png" width="50%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/33.png" width="45%">  
  
For the vector **b** we will represent the right side **f** as an interpolating function:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/38.png" width="35%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/36.png" width="45%">  

### Global Matrix
Now, using the local matrices and the local vector, we need to assemble the global matrix (**A'**) and the global vector of the right side (**f**). Consider the following area as an example:   
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/50.png" width="60%">  
This picture shows the finite element numbers and global nodes. To find the solution we need to calculate 5 local matrices and 5 local vectors. The global matrix will be assembled from the local matrices: the assembly takes place in accordance with the local node numbering on the finite elements. For the second finite element, for example, local number **1** corresponds to global number **3**, local number **2** corresponds to global number **4**, and so on. Let's add the first two local matrices to the global matrix:    
  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/47.png" width="75%">  
  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/49.png" width="75%">  
  
The global vector can be found in the same way.    

### Sparse String Format
Since each of the elements is related to a limited number of neighboring elements, the system of linear algebraic equations has a [sparse form](https://en.wikipedia.org/wiki/Sparse_matrix). We will store the global matrix efficiently, in the sparse string format:  
- **ggl**, **ggu** - arrays for non-zero elements of the lower and upper triangles
- **diOrig** - array for diagonal elements
- **ig** - array for indexes of elements from which a string of non-zero elements begins
- **jg** - array for the column/row numbers in which the lower/upper triangle element is located, respectively
  
This allows us to store only non-zero elements. Consider the matrix:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/51.png" width="55%">  
Its representation in the sparse string format:
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/52.png" width="70%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/53.png" width="33%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/54.png" width="26%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/55.png" width="27%">  

### Boundary Conditions 
Up to this point, boundary conditions have not been taken into account, so it's time to do it.  
  
To account for natural boundary conditions, we form local matrices and local vectors that will be added to the linear system similarly to the local matrices and local vectors of finite elements. In our (two-dimensional) case, the natural boundary conditions are considered on edges (one-dimensional segments).  
Let's denote by **Г** the edge of length **h**, and by **(i, j)** the edge coordinates. We will represent the parameters **θ** and **u** as linear basis expansions, the parameter **β** will be regarded as a constant. Only two basis functions (**Ψ1** and **Ψ2**) are nonzero on each edge, so the corresponding expressions take the following form.  
  
For boundary condition of the second type:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/78.png" width="40%">  
For boundary condition of the third type:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/80.png" width="36%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/81.png" width="49%">  
  
Taking the essential (first type) boundary conditions into account, which is carried out after the global matrix is assembled, is the follows. By going through all the nodes on the edge, we replace the corresponding elements of the diagonal of the global matrix with a number that is much larger than the rest of the matrix elements.  
In the global vector, the element with the corresponding number is assigned the product of the large number and the analytic value of the function.  

### LU Decomposition and Gaussian Elimination
We will decompose the global matrix into lower and upper triangles:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/56.png" width="65%">  
using the formulas:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/57.png" width="49%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/58.png" width="52%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/59.png" width="29%">  

[Lower-upper decomposition](https://en.wikipedia.org/wiki/LU_decomposition) (or factorization) is a better way to implement [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination), especially for repeated solving a number of equations with the same left-hand side. The original system can now be solved in two steps:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/70.png" width="25%">  

### Local Optimal Scheme
To solve the final system of linear equations, we will use a method called the Local Optimal Scheme (LOS) with [incomplete factorization](https://en.wikipedia.org/wiki/Incomplete_LU_factorization), which is a lesser-known alternative to the [conjugate gradient method](https://en.wikipedia.org/wiki/Conjugate_gradient_method) (CGM).   
  
Incomplete factorization allows to solve the following equivalent system instead of solving **LUq=f**:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/60.png" width="25%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/61.png" width="20%">  
Now let's assume:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/62.png" width="50%">  
Then we start the iterative process:  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/63.png" width="30%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/64.png" width="30%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/65.png" width="30%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/66.png" width="40%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/67.png" width="34%">  
<img src="https://github.com/lenferdetroud/misc/blob/master/fem/68.png" width="39%">  
The process ends when <img src="https://github.com/lenferdetroud/misc/blob/master/fem/69.png" width="7%"> becomes small enough. Vector **x** is the solution (i.e., it contains the coefficients **q1**, **q2**, ..., **qn**).  

## Implementation
- [main.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/main.cpp)  
The main program module containing function calls, local matrices and local vectors construction, diffusion coefficient decomposition, global matrix assembling and linear system solution.  
- [config.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/config.cpp)  
A module designed to configure the parameters of the area, the source function, the vector of the right side and the LOS parameters. This module also allocates memory, configures pointers, and builds the matrix portrait.  
- [boundaries.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/boundaries.cpp)  
The implementation of a sequential accounting of natural and essential boundary conditions by enumerating the given types of conditions on the boundaries of the domain.  
- [gauss.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/gauss.cpp)   
The implementation of forward and backward Gaussian elimination, as well as global matrix decomposition. These functions allow the implementation of incomplete factorization for the LOS.
- [math.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/math.cpp)  
A module containing minor functions: multiplication of global matrix by vector and vector by vector (applied in LOS), sum of products of upper and lower triangles (applied in global matrix decomposition). 
- [io.cpp](https://github.com/lenferdetroud/finite-element-method/blob/main/io.cpp)    
A module containing functions to read input data and build a table with the result of program execution on the given data.
- [common.h](https://github.com/lenferdetroud/finite-element-method/blob/main/common.h)  
A header file containing preprocessor directives, function declarations and global variables common to all modules.
- [size.txt](https://github.com/lenferdetroud/finite-element-method/blob/main/size.txt)   
This file specifies the number of nodes on the x-axis and the number of nodes on the y-axis, separated by a space.
- [nodesX.txt](https://github.com/lenferdetroud/finite-element-method/blob/main/nodesX.txt), [nodesY.txt](https://github.com/lenferdetroud/finite-element-method/blob/main/nodesY.txt)  
Files containing node coordinates in the X and Y axes, respectively.
- [edgeTypes.txt](https://github.com/lenferdetroud/finite-element-method/blob/main/edgeTypes.txt)   
A file containing four space-separated values (1/2/3) describing the type of condition on the boundaries of an area. The first value corresponds to the lower boundary, the second to the right, the third to the top, and the fourth to the left.

## Testing and Research
You can find all the tests in the file [tests.pdf](https://github.com/lenferdetroud/finite-element-method/blob/main/tests.pdf).

## Install
Since this is a Windows console application, you can simply clone the repository directly into your Visual Studio environment or use `git clone https://github.com/lenferdetroud/finite-element-method.git`. The application has no interface.  
