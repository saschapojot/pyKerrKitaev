import numpy as np

#fermion number
N=300
kIndAll=range(0,N)
dk=2*np.pi/N

#before quench
mu0=3
t0=1
d0=1

#after quench
mu1=3
t1=t0
d1=d0

#nonlinearity strength
lmd=0

#small time step number
R=20
#large time step number
Q=100

#small time step
dt=0.01
#large time step
ds=R*dt

tol=1e-16
cutOff=1.2

threadNum=12
