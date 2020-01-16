'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given an Unsafe Set and Matrix A, compute the amount of error that can be
introduced in A and still be safe.

Documentation: Not yet available. (TODO)
'''

from BloatAPI import *
import math
import scipy.special as SP
from mpmath import *

class RobustMetric:
    '''
    Computes Robustness Metric according to the following
    three Bloating Methods:
        1. Kagstrom 1
        2. Kagstrom 2
        3. Loan
    '''

    def __init__(self,A,U,m,IS,t):
        self.A=A # The Linear Dynamical Matrix A
        self.unsafe=U # The Unsafe Set
        '''
        U is a numpy array of intervals
        U[i,0]=(a,b)
        '''
        self.method=m # Method to be used to compute Robustness Metric (Kagstrom1, Kagstrom2, Loan)
        self.initialSet=IS # Initial Set
        '''
        initialSet[i,j]=mpi(a,b)
        '''
        self.t=t # Time

    def getSafeBloat(self):
        reachSet=np.matmul(SLA.expm(self.A*self.t),self.initialSet)
        min=9e100
        ind=-1
        for i in range(self.A.shape[0]):
            a1=float(nstr(reachSet[i][0]).split(',')[0][1:])
            b1=float(nstr(reachSet[i][0]).split(',')[1][:-1])
            a2=self.unsafe[i][0][0]
            b2=self.unsafe[i][0][1]
            if a2>=a1 and a2<=b1:
                print("Already Unsafe!!")
                return -9999
            elif a1>=a2 and a1<=b2:
                print("Already Unsafe!!")
                return -9999
            elif b1<=a2:
                dis=abs(a2-b1)
            elif b2<=a1:
                dis=abs(a1-b2)
            if dis<min:
                min=dis
                ind=i
            dis=0
        a1=float(nstr(reachSet[ind][0]).split(',')[0][1:])
        b1=float(nstr(reachSet[ind][0]).split(',')[1][:-1])
        a2=self.unsafe[ind][0][0]
        b2=self.unsafe[ind][0][1]
        c=abs(b1-a1)/2
        rad=b1-c
        #print(rad,c,ind)
        if a2>b1:
            return abs(a2-c)/rad
        elif b2<a1:
            return abs(b2-c)/rad

    def getRobustMetric(self):
        if self.method.lower()=='kagstrom1':
            return self.getRobustMetricKagstrom1()
        elif self.method.lower()=='kagstrom2':
            return self.getRobustMetricKagstrom2()
        elif self.method.lower()=='loan':
            return self.getRobustMetricLoan()

    def getRobustMetricKagstrom1(self):
        beta=self.getSafeBloat()
        kag1=BloatKagstrom(self.A,[])
        (Q,lam,M)=kag1.decompose()
        normM=IntervalNorm.spectralNorm(M)
        Mfac=BloatKagstrom.computeP(normM*self.t,self.A.shape[0])
        robMet=math.log((beta/Mfac)+1)/(Mfac*self.t)
        return robMet


    def getRobustMetricKagstrom2(self):
        beta=self.getSafeBloat()
        kag2=BloatKagstrom(self.A,[])
        (S,N,l,ep)=BloatKagstrom.JNFDecomp(self.A)
        D=BloatKagstrom.getD(l,ep,N)
        KSD=BloatKagstrom.computeK(np.matmul(S,D))
        e_ep=np.exp(ep*self.t)
        robMet=math.log((beta/(KSD*e_ep))+1)/KSD*self.t
        return robMet

    def getRobustMetricLoan(self):
        beta=self.getSafeBloat()
        loan=BloatLoan(self.A,{})
        e_pow=np.exp(math.log(beta)-(loan.computeMu()*self.t)+(loan.computeAlpha()*self.t))
        robMet=SP.lambertw(e_pow.real)/self.t
        return robMet.real

    def getMaxSafeBloat(self):
        print("Started...")
        b1=self.getRobustMetricKagstrom1()
        b2=self.getRobustMetricKagstrom2()
        b3=self.getRobustMetricLoan()
        print("...Done")
        b=max(b1,b2,b3)
        if (b1==b):
            return (b,'Kagstrom1')
        elif (b2==b):
            return (b,'Kagstrom2')
        else:
            return (b,'Loan')

# Tester
if True:
    A=np.array([
    [0,0,1,0,1],
    [1,0,0,1,1],
    [0,0,1,0,0],
    [0,0,1,0,1],
    [0,1,1,0,0],
    ])
    U=np.array([
    [(-20000000000000000000000000,-10000000000000000000000000)],
    [(-20000000000000000000000000,-10000000000000000000000000)],
    [(-20000000000000000000000000,-10000000000000000000000000)],
    [(-20000000000000000000000000,-10000000000000000000000000)],
    [(-20000000000000000000000000,-10000000000000000000000000)],
    ])
    '''U=np.array([
    [(-7,-6)],
    [(-7,-6)],
    [(-7,-6)],
    [(-7,-6)],
    [(-7,-6)],
    ])'''
    IS=np.array([
    [mpi(1,3)],
    [mpi(1,3)],
    [mpi(1,3)],
    [mpi(1,3)],
    [mpi(1,3)],
    ])
    t=20
    method='Loan'

    rM=RobustMetric(A,U,method,IS,t)
    print(rM.getSafeBloat())
    #print(rM.getRobustMetric())
    print(rM.getMaxSafeBloat())
