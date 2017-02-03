import numpy as np
import matplotlib.pyplot as plt

# Actual computation
nx = 20; nt = 50
dt = 0.01; c = 1
dx = 2/(nx - 1)
nu = 1
import numpy as np
x = np.linspace(0,2,nx)
u = np.zeros(nx)

for i in range(nx):
    if x[i]>=0.5 and x[i]<=1:
       u[i] = 2
    else:
       u[i] = 1

for it in range(1,nt+1):
    un = u.copy()
    for i in range(1,nx-1):
       u[i] = un[i] + nu*(dt/dx**2)*(un[i+1] - 2*un[i] + un[i-1])

print("After "+ str(nt)+ " time steps,the solution is :\n")
print(u)
import matplotlib.pyplot as plt
plt.plot(x,u, 'b')
plt.grid(True)
plt.xlabel("x")
plt.ylabel("U")
plt.title("1D diffusion equation solution after 50 time steps")
plt.show()
