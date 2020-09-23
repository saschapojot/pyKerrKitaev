
import scipy.linalg as slin
from datetime import datetime
from multiprocessing import Pool
import matplotlib.pyplot as plt
from  initDat import *

def b0(k):
    return d0 * np.sin(dk * k)


def c0(k):
    return t0 * np.cos(dk * k) + mu0 / 2

#initial vector with kth momentum
def initVec(k):
    b0Val = b0(k)
    c0Val = c0(k)
    denom2 = 2 * b0Val ** 2 + 2 * c0Val ** 2 - 2 * c0Val * np.sqrt(b0Val ** 2 + c0Val ** 2)
    if denom2 > tol:
        denom = np.sqrt(denom2)
        return np.array([1j * b0Val / denom, (c0Val - np.sqrt(b0Val ** 2 + c0Val ** 2)) / denom])
    else:
        if c0Val>0:
            return np.array([1,0])
        else:
            return np.array([0,1])

def H0(k):
    '''

    :param k:
    :return: linear part of Hamiltonian after quench
    '''
    tmp00=t1*np.cos(k*dk)
    rst=np.zeros((2,2),dtype=complex)

    rst[0,0]=-mu1-tmp00
    tmp01=d1*np.sin(k*dk)
    rst[0,1]=1j*tmp01
    rst[1,0]=-1j*tmp01
    rst[1,1]=tmp00
    return rst
#calculate each H0(k)
for k in kIndAll:
    H0All.append(H0(k))
def S2(k, vecStart):
    '''
    One step S2
    :param k:
    :param vecStart:
    :return:
    '''
    h0Val=H0All[k]
    exph0Val=np.matrix(slin.expm(-1j/2*dt*h0Val))

    #step 1
    vec1=np.array(exph0Val.dot(vecStart))
    #step 2

    vkm=vec1[0][0]
    wkm=vec1[0][1]

    etakm=vkm*np.exp(-1j*lmd*np.abs(vkm)**2*dt)
    zetakm=wkm*np.exp(-1j*lmd*np.abs(wkm)**2*dt)
    #step 3
    vec2=np.array([etakm,zetakm])
    vec3=np.array(exph0Val.dot(vec2))

    return np.array([vec3[0][0],vec3[0][1]])

def calculateVec(k):
    '''
    calculate the state vectors starting with the kth momentum
    :param k:
    :return:
    '''
    #initialization
    veck0=initVec(k)
    statesFromK=[]
    statesFromK.append(veck0)
    #for time step m=0,1,...,Q*R-1
    for m in range(0,Q*R):
        vecCurr=statesFromK[-1]
        vecNext=S2(k,vecCurr)
        statesFromK.append(vecNext)

    return [k,statesFromK]




pool1=Pool(threadNum)
retAll=pool1.map(calculateVec,kIndAll)
pool1.close()
pool1.join()
for item in retAll:
    kVal=item[0]
    statesAll[kVal]=item[1]


def Jkab(k,a,b):
    '''

    :param k: 0,1,...,N-1
    :param a: 0,1,...,Q-1
    :param b: 0,1,...,R
    :return:
    '''
    vec=statesAll[k][a*R+b]
    HVal=np.zeros((2,2),dtype=complex)
    yVal=vec[0]
    zVal=vec[1]
    tmp00=t1*np.cos(k*dk)
    tmp01=d1*np.sin(k*dk)
    HVal[0,0]=-mu1-tmp00+lmd*np.abs(yVal)**2
    HVal[0,1]=1j*tmp01
    HVal[1,0]=-1j*tmp01
    HVal[1,1]=tmp00+lmd*np.abs(zVal)**2
    tmp=np.conj(vec.transpose()).dot(HVal).dot(vec)
    tmp/=(np.abs(yVal)**2+np.abs(zVal)**2)

    return tmp

def simpsonD(k,a):
    '''

    :param k: 0,1,...,N-1
    :param a: 0,1,...,Q-1
    :return: integral from ads to (a+1)ds
    '''
    oddVals=[]
    evenVals=[]
    for b in range(1,R,2):
        oddVals.append(Jkab(k,a,b))
    for b in range(2,R,2):
        evenVals.append(Jkab(k,a,b))

    tmp=Jkab(k,a,0)+4*sum(oddVals)+2*sum(evenVals)+Jkab(k,a,R)
    tmp*=dt/3
    return tmp

def oneSimpTabEntry(kAndAOnePair):
    kVal=kAndAOnePair[0]
    aVal=kAndAOnePair[1]
    return [kVal,aVal,simpsonD(kVal,aVal)]

kAndAAllPairs=[[kVal, aVal] for kVal in kIndAll for aVal in range(0,Q)]
pool2=Pool(threadNum)
KASimpD=pool2.map(oneSimpTabEntry,kAndAAllPairs)
pool2.close()
pool2.join()
for item in KASimpD:
    kVal=item[0]
    aVal=item[1]
    tmp=item[2]
    simpTab[kVal][aVal]=tmp

def thetaDTabOneEntry(kAndQPair):
    '''

    :param k: 0,1,...,N-1
    :param q: 1,2,...,Q
    :return:
    '''
    sumTmp=0
    kVal=kAndQPair[0]
    qVal=kAndQPair[1]
    for a in range(0,qVal):
        sumTmp+=simpTab[kVal][a]
    sumTmp*=-1

    vec0=statesAll[kVal][0]
    vecCurr=statesAll[kVal][qVal*R]

    denomTmp=np.conj(vec0.transpose()).dot(vec0)
    numerTmp=np.conj(vecCurr.transpose()).dot(vecCurr)
    return [kVal,qVal,sumTmp+0.5*1j*np.log(numerTmp/denomTmp)]

#for k=0,1,...,N-1 q=0, thetaD=0
for kVal in kIndAll:
    thetaDTab[kVal][0]=0
#for k=0,1,...,N-1, q=1,2,...,Q
kAndQAllPairs=[[kVal, qVal] for kVal in kIndAll for qVal in range(1,Q+1)]

pool3=Pool(threadNum)
KQDVals=pool3.map(thetaDTabOneEntry,kAndQAllPairs)
pool3.close()
pool3.join()
for item in KQDVals:
    kVal=item[0]
    qVal=item[1]
    tDVal=item[2]
    thetaDTab[kVal][qVal]=tDVal


def thetaTotTabOneEntry(kAndQOnePair):
    '''

    :param kAndQOnePair: k=0,1,...,N-1;q=0,1,...,Q
    :return:
    '''
    kVal=kAndQOnePair[0]
    qVal=kAndQOnePair[1]

    vecCurr=statesAll[kVal][qVal*R]
    vec0=statesAll[kVal][0]

    numerTmp=np.conj(vec0.transpose()).dot(vecCurr)

    tmp=-1j*np.log(numerTmp/np.abs(numerTmp))
    return [kVal,qVal,tmp]

kAndQForThetaTotAllPairs=[[kVal, qVal] for kVal in kIndAll for qVal in range(0,Q+1)]
pool4=Pool(threadNum)
KQToTAll=pool4.map(thetaTotTabOneEntry,kAndQForThetaTotAllPairs)
pool4.close()
pool4.join()
for item in KQToTAll:
    kVal=item[0]
    qVal=item[1]
    tTotVal=item[2]
    thetaTotTab[kVal][qVal]=tTotVal

#fill thetaGTab
for k in kIndAll:
    for q in range(0,Q+1):
        thetaGTab[k][q]=thetaTotTab[k][q]-thetaDTab[k][q]

def jumpDecision(incr):
    tmp=incr/np.pi
    if tmp>=cutOff:
        return incr-2*np.pi
    elif tmp<=-cutOff:
        return incr+2*np.pi

    else:
        return incr
#fill beta
for q in range(0,Q+1):
    for k in range(0,N-1):
        incrTmp=np.real(thetaGTab[k+1][q]-thetaGTab[k][q])
        beta[q][k]=jumpDecision(incrTmp)

#fill W
for q in range(0,Q+1):
    tmp=0.0
    for k in range(0,int(N/2)):
        tmp+=beta[q][k]
    tmp/=2*np.pi
    W.append(tmp)



Ts=np.arange(0,Q+1)*ds
plt.figure()
plt.plot(Ts,W)
plt.yticks(np.arange(min(W),max(W)+1))
plt.show()

plt.close()