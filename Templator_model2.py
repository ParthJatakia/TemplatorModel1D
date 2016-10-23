import numpy as np
import matplotlib.pyplot as plt


r = 0.775
ku = 0.015874
Dv = 6.29961
Du = 31.498
K = 0.15874

size = 128  # size of the 2D grid
dx = 2. / size  # space step
T = 1.0  # total time
dt = dx**2/2  # time step
n = 300

U = [2]*size
V = [None]*0
for i in range(size):
    if i > 40 and i < 85:
        V.append(np.sin(i)+1)
    else:
        V.append(0)
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
    print "reached %i" %i
    # We take the values of u and v inside the grid.
    Uc = U[1:-1]
    Vc = V[1:-1]
    # We update the variables.
    print len(U)
    print len(V)
    U[1:-1], V[1:-1] = \
        Uc + dt * (r - ku*np.square(Uc) -np.multiply(np.square(Uc),Vc)+ Du*deltaU), \
        Vc + dt * (ku*np.square(Uc)+ np.multiply(np.square(Uc), Vc) -np.divide(Vc, np.add(K,Vc))+Dv*deltaV)
    print len(U)
    print len(V)
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
plt.xticks([]); plt.yticks([])
plt.colorbar()
plt.show()
