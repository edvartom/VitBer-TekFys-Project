The first figure shows the maximum error as a function of tolerance. From what we see in the plot, it seems that the error increases linearly for larger tolerance. This is expected, as a larger tolerance implies that the simulation will tolerate larger absolute differences between the calculated values of $y$ for our second order ode, and for the third order ode.

<br>

The second figure shows the number of timesteps used by Runga-Kutta function, uses given some value $\alpha$, where $\alpha$ is the pessimistic factor. This factor ensures that $h_{new}$ is smaller than the current $h$, ensuring that the calculated value $y$ falls within the tolerance. In the figure, we see that the gradient is steep at first, but decreases quickly as $\alpha$ is increased, meaning that the number of timesteps required falls off quickly for small $\alpha$ as $\alpha$ is increased. This is because $\alpha \propto h$, so that when alpha is increased, the size of the timestep is increased. As such, the number of timestep needed falls.

Almost all of the timesteps are accepted for small values of $\alpha$. For larger $\alpha$, fewer steps are accepted, e.g. rejected. It therefore makes sense that the trend flattens out (and does not go down linearly) as $\alpha$ increases. 

As we can see in the figure, it reaches a minimum around $\alpha = 0.9$. The minimum occurs when the stepsize large enough to only require a small number of them are needed, while also being small enough to accept the majority of the timesteps. 

For $\alpha$ close to $1$, the gradient is positive because the timesteps sizes are too large. Steps must therefore be discarded, and recalculated, resulting in the Runga Kutta requiring additional steps. 

In conclusion, to optimise the simulation, we want to minimise the number of steps, and tehrefore, should choose $\alpha \approx 0.9$. 