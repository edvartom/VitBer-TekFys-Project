from funcs import *

def g(x):
    return x + np.sin(x) + np.cos(x)

x = np.linspace(-5,5,1000)
z = secant_method(-1, 4, g)

print("g(z) =", g(z))
print("z =", z)

plt.plot(x, g(x), label="g(x)")
plt.plot(z, g(z), 'ro', label="Root ="+str(round(z,3)))
plt.title("Root of $g(x)=x+sin(x)+cos(x)$")
plt.ylabel("$g(x)$")
plt.xlabel("$x$")
plt.grid()
plt.legend()
plt.show()