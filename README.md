# Partial Differential Equations 3 (SCEE09004)
This repository contains Jupyter workbooks for the Partial Differential Equations 3 course.  The workbooks cover the **numerical methods** half of the course and use finite difference methods for the solution of eliptic, parabolic and hyperbolic differential equations.  All of these workbooks make use of the <code>matplotlib</code>, <code>numpy</code> and <code>scipy</code> python libraries.

Each workbook covers a different topic from the course.  Many of them make use of the <code>refinement_analysis.py</code> module which must be in the same folder as the workbook so it can be found by the python import:

<code>from refinement_analysis import refinement_analysis</code>

These workbooks contain *Worked examples*.  These are complete solutions with explanations of work being conducted.  Separate *laboratory problems* will be issued which are incomplete and will involve both mathematical and numerical methods for the solution of a particular problem. These will contain the fundamental python tools needed to contruct a numerical solution but will need Python code to be written.

All code is covered by the Creative-Commons by Attribution Licenses and is copyright to The School of Engineering, University of Edinburgh.  The code was written by Professor David Ingram.

## Worked examples
1. Laplace Itterative - Jacobi, Gauss-Siedel and SOR methods for the solution of the Laplace equations.
2. Laplace Faster - Making the Gauss-Siedel and SOR methods more computationally efficient by using checker-boarding and other techniques to speed up loops.
3. V&V - Verification and validation of solvers.  Testiting and checking accuracy including mesh refinement studies (needs <code>refinement_analysis.py</code>).
4. Boundary Conditions - Implementation of Dirichlet, Neumann and Robin boundary conditions for the Laplacian (needs <code>refinement_analysis.py</code>).
5. BiCGStab - The Bi-stabilised conjugate gradient method: A fast, itterative, implicit matrix solver to be used in place of the itterative methods we've used previously (needs <code>refinement_analysis.py</code>).
6. Source Terms - Solving the Poisson equation [i.e. the inhomogenous version of the Laplacian] (needs <code>refinement_analysis.py</code>). 
7. Explicit Parabolic Solver - The Forward Time-Centered Space (FTCS) solver for the solution of 1-D parabolic equations (needs <code>refinement_analysis.py</code>).
9. Implicit Parabolic Solver - Using the Crank-Nicholson method to solver for the solution of 1-D parabolic equations (needs <code>refinement_analysis.py</code>).
10. 2D Parabolic - Solving 2D parabolic problems using opperator splitting (the solver used is the explicit FTCS scheme, though it could be modified to use Crank-Nicholson).
11. Hyperbolic Solver - Using the Lax-Wendroff and First Order Upwind schemes to solve the 1D linear advection problem with periodic boundary conditions.  

## Laboratories
The Laboratories folder includes a Jupyter workbook for each of the Laboratories from weeks 1 to 10.

## Coursework
The couursework folder contains Jupyter notebooks needed for the coursework assignments.
