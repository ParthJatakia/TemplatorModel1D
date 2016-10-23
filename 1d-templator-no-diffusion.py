import numpy as np
import matplotlib.pyplot as plt


r = 0.523
ku = 0.01
Dv = 0
Du = 0
K = 0.3

dt = 0.1  # time step
n = 5000

U = [0.0]*n
V = [0.0]*n

U[0]=1
V[0]=1

# We simulate the PDE with the finite difference method.

for i in range(n-1):
    # We take the values of u and v inside the grid.
    Uc = U[i]
    Vc = V[i]
    # We update the variables.
    U[i+1] = Uc + dt*( r - ku*(Uc**2) - (Uc**2)*Vc )
    V[i+1] = Vc + dt*( (ku*(Uc**2) + (Uc**2)*Vc) - (Vc/(K+Vc)))


plt.plot(range(n),V)
plt.show()