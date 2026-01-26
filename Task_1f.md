## 1f) 

In this task, we will solve the boundary value problem

$$
\frac{d^2y}{dx^2} = -4 \cdot \sin(2x), \quad y(0) = 0, \quad y(2\pi) = 0.
$$

Just like before, we define 

$$
\vec{y} :=
\begin{bmatrix}
y(x) \\
y'(x)\\
\end{bmatrix},
$$

and 
$$
f(\vec{y}, x) := \frac{d\vec{y}}{dx} = \begin{bmatrix}
y'(x) \\
-4\cdot\sin(2x)\\
\end{bmatrix}
.$$

We must now find a $b:=y'(0)$ so the following boundary conditions 
$$
y_b(0)=y_b(2\pi)=0 \tag{1}
$$ 

is satisfied. 

The secant method finds the root, e.g. the $x$-value for which $y(x)=0$. Given that the the boundary condition at the end for $y(x)$ is $0$, we may use the Runga-Kutta method to approximate a solution for $y(x)$. 

The process is as follows. We initiate a guess for $b$, and use the Runga-Kutta method to approximate to the IVP condition $y'(x) = b$, which returns $y(x_{end})$. The secant method will repeat in a loop till some value of $b$ that satisifes (1) is found, within a given tolerance. 


# Placeholder, here is the code

The figure shows the solutions to the BVP. The BVP is solved as an IVP using the Runga-Kutta mehtod for varying values of $b$, where $b$ is found using the secant_method function. The low tolerance in $b$ leads to the last two $b$-values being similar enough that their plots overlap in the figure. 