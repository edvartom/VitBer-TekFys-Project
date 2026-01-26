## Exercise 1
### 1a)

We want to solve the initial value problem 

\frac{d^2y}{dx^2} = -4\cdot \sin(2x), \quad y(0)=0, \quad y'(0) = 2\text{.}$$

We begin by integrating both sides twice. First, we obtain

$$
\frac{dy}{dx} = 4\cdot \frac{1}{2} \cos(2x) + C_1 =2\cos(2x) + C_1 \text{.}
$$

Thereafter, we obtain

$$
y(x) = \sin(2x) + C_1x + C_2 \text{.}
$$

The first initial condition gives us 

$$
y(0) = \sin(2\cdot 0) + C_1\cdot 0 + C_2 \Rightarrow C_2 = 0.
$$

We know the derivative from above, and put in the second initial condition as well. We get 

$$
y'(0) = 2\cos(2\cdot 0) + C_1 = 2 \Rightarrow 2\cdot 1 + C_2 = 2 \Rightarrow C_2 = 0 \text{.}
$$

This gives the solution

$$
y(x) = \sin(2x).
$$
