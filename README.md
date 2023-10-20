# Solution of the initial boundary value problem for the heat equation

## The initial boundary value problem for a nonlinear one-dimensional heat equation with a source is considered:
```math
\partial u/ \partial t = \partial / \partial x (k(u) \partial u / \partial x) + f(t,x,u),  x \in (0,1), t > 0; \  

u(0, x) = \phi(x), x \in [0,1];    
\alpha_0 u(t,0) + \beta_0u_x(t,0) = \psi_0(t), t > 0    
 \alpha_1 u(t,1) + \beta_1u_x(t,1) = \psi_1(t), t > 0    
```
The linear case of the heat equation is considered
```math
k(u)=k_0=Const, f(t,x,u)=f(t,x).  
```
Next problems are solved in this program:
Problem 1. Explicit difference scheme. (left corner)
Problem 2. Implicit difference scheme. (Krank-Nicholson scheme)
Problem 3. Solving a nonlinear problem using a conservative scheme on a uniform grid.
