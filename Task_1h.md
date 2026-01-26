# 1h)

In the first figure, we have plotted the solution to eq. (2) for both the Scipy BVP solver module,  and the IVP solver and root finder from 1g). The solutions appears to be identical. However, the Scipy BVP module was faster than our solver by 12.6s. 

The second figure shows the absolute differences between the two BVP solvers as a function of $x$. The The difference appears to be of order $10^{-3}$ order smaller than the actual values of $y(x)$, and is therefore neglectible. 

In conclusion, by considering the time difference between Scipy's BVP solver and our solver, it is optimal to use Scipy's BVP solver. While the time difference is not considerably large in this problem, we should expect a longer runtime for more complicated problems that require more iterations. 