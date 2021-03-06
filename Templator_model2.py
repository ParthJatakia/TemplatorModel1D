import numpy as np
import matplotlib.pyplot as plt


r = 0.775#0.775
ku = 0.015874
Dv = 0.015/5
Du = 0.075/5
K = 0.15874

size = 200  # size of the 1D grid
dx = 2. / size  # space step
l=size/2
#T = 1.0  # total time
dt = 0.001 #dx**2/100  # time step
n = 50000

U = [1]*size
V = [None]*0
for i in range(size):
    if i > (size - l) / 2 and i < (size + l) / 2:
        V.append(np.cos(3 * np.pi * (i - (size - l) / 2) / l) ** 2)
    else:
        V.append(1.0)
Ut = U

Vt = V


def laplacian(Z):
    # Ztop = Z[0:-2,1:-1]
    Zleft = Z[0:-2]
    # Zbottom = Z[2:,1:-1]
    Zright = Z[2:]
    Zcenter = Z[1:-1]
    return np.divide(np.subtract(np.add(Zleft,Zright), np.multiply(Zcenter,2)),dx**2)

# We simulate the PDE with the finite difference method.

for i in range(n):
    # We compute the Laplacian of u and v.
    deltaU = laplacian(U)
    deltaV = laplacian(V)
    # We take the values of u and v inside the grid.
    Uc = U[1:-1]
    Vc = V[1:-1]
    # We update the variables.
    U[1:-1], V[1:-1] = \
        Uc + dt * (r - ku*np.square(Uc) - np.multiply(np.square(Uc),Vc)+ Du*deltaU), \
        Vc + dt * (ku*np.square(Uc)+ np.multiply(np.square(Uc), Vc) -np.divide(Vc, np.add(K,Vc))+Dv*deltaV)
    if i%(n/50)==0:
        Ut = np.vstack((U,Ut))
        Vt = np.vstack((Vt, V))


    # Neumann conditions: derivatives at the edges
    # are null.
    for Z in (U, V):
        Z[0] = Z[1]
        Z[-1] = Z[-2]
        #Z[:,0] = Z[:,1]
        #Z[:,-1] = Z[:,-2]


plt.imshow(Vt, extent=[-1,1,-1,1], cmap='spectral')
plt.xlabel('position')
plt.ylabel('time (increases downwards)')
plt.colorbar()
plt.show()
