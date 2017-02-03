# 1D Burgers equation
# du/dt + u*du/dx = v*d^2u/dx^2
# Discretisation using FD in t, BD for convective, CD for diffusion

from numpy import *
from matplotlib.pyplot import *

# Actual computation
nx = 20; nt = 100
dt = 0.01
dx = 2/(nx - 1)
nu = 0.3
import numpy as np
x = linspace(0,2,nx)
u = zeros(nx)

for i in range(nx):
    if x[i]>=0.5 and x[i]<=1:
       u[i] = 2
    else:
       u[i] = 1
uzero = u.copy()

line, = plot(x,uzero, 'b--')
grid(True)
#axis([0, 2, 0, 2])
xlabel("x")
ylabel("U")
title("1D burgers equation solution after "+str(nt)+" time steps")
for it in range(1,nt+1):
    un = u.copy()
    for i in range(1,nx-1):
       u[i] = un[i] - un[i]*(dt/dx)*(un[i] - un[i-1]) + nu*(dt/dx**2)*(un[i+1] - 2*un[i] + un[i-1])

    line.set_ydata(u)
show()
hold(True)  
# Plotting
print("After "+ str(nt)+ " time steps,the solution is :\n")
print(u)

