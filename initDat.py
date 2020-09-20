#dict to hold all state vectors
#statesAll[k]=[vec0,vec1,...,vec_{QR}]
from consts import *
statesAll=dict()
H0All=[]


#simpTab is an array,k=0,1,...,N-1
#a=0,1,...,Q-1
simpTab=[]
for k in range(0,N):
    simpTab.append(list(np.arange(0,Q)))
#thetaDTab is an array,k=0,1,...,N-1
#q=0,1,...,Q
thetaDTab=[]
for k in range(0,N):
    thetaDTab.append(list(np.arange(0,Q+1)))


#thetaTotTab is an array, k=0,1,...,N-1
#q=0,1,...,Q
thetaTotTab=[]
for k in range(0,N):
    thetaTotTab.append(list(np.arange(0,Q+1)))

#thetaGTab is an array, k=0,1,...,N-1
#q=0,1,...,Q
thetaGTab=[]
for k in range(0,N):
    thetaGTab.append(list(np.arange(0,Q+1)))

#beta is an array, q=0,1,...,Q
#k=0,1,...,N-2
beta=[]
for q in range(0,Q+1):
    beta.append(list(np.arange(0,N-1)))

W=[]