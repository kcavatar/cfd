# 2D convection equation
# du/dt + u*du/dx + v*du/dy = 0
# dv/dt + u*dv/dx + v*dv/dy = 0


import numpy as np
import matplotlib.pyplot as plt

# Actual computation
nx = 20; ny = 20
nt = 50
dt = 0.01; c = 1
dx = 2/(nx - 1); dy = 2/(ny - 1)

import numpy as np
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
u = np.zeros(nx)
v = np.zeros(ny)
for i in range(nx):
    if (x[i]>=0.5 and x[i]<=1) and (y[i]>=0.5 and y[i]<=1):
       u[i] = 2
       v[i] = 2
    else:
       u[i] = 1
       v[i] = 1

for it in range(1,nt+1):
    un = u.copy()
    vn = v.copy()
    for i in range(1,nx-1):
       u[i] = un[i] - un[i]*(dt/dx)*(un[i] - un[i-1]) - vn[i]*(dt/dy)*(un[i] - un[i-1])
       v[i] = vn[i] - un[i]*(dt/dx)*(vn[i] - vn[i-1]) - vn[i]*(dt/dy)*(vn[i] - vn[i-1])
# Results
print("After "+str(nt)+" time steps, the values are as :\n")
print("x - component of velocity, u :")
print(u)
print("y - component of velocity, v :")
print(v)
