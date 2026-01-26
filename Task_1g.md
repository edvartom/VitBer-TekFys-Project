# Suggested Reformulation

## 1g)
We proceed by to solve the BVP given by 
$$
\frac{d^2y}{dx^2} = y(x)+\sin(x), 
\quad y(0)=y(12)=0 \tag2
$$

using the IVP-solver and root finder from the previous tasks.  

First, we repeat the steps in 1f); define a 2D vector, $\vec{y}$ given by 

$$
\vec{f}=
\begin{bmatrix}
y'(x) \\
y(x) + sin(x)
\end{bmatrix} \text{.}
$$

Hence, reformulate our BVP into IVP that may solved numerically. Our boundary conditions are 

$$
y_b(0)=y_b(12)=0 \text{.}
$$

<br>


# Placeholder, put code block here

The graph shows the numerically calulated solution for y(x) plotted for x, given the BVP described in (2). The solution appears to be an altered negative sinusodial wave, which is also what you would get when integrating (2) twice. # TODO: Dobbeltsjekk dette