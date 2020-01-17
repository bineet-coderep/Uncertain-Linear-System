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

    def __init__(self,A,B,U,m,IS,t,mode,h):
        self.A=RobustMetric.createMatrix(A,B,mode,h) # The Linear Dynamical Matrix A
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

    @staticmethod
    def createMatrix(A,B,mode,h=0):
        ''' Creates a single matrix based on
        . or +.
        In case of . a roungh approximation is
        done'''

        n1=np.size(A,0)
        if (np.size(B)>0):
            n2=np.size(B,1)
        else:
            n2=0
        n=n1+n2

        C=np.zeros((n,n))
        for i in range(n1):
            for j in range(n1):
                C[i][j]=A[i][j]
        for i in range(n1):
            j2=0
            for j in range(n1,n1+n2):
                C[i][j]=B[i][j2]
                j2=j2+1

        if mode=='+':
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C

    def getSafeBloat(self):
        min=9e100
        ind=-1
        reachSet=np.matmul(SLA.expm(self.A*self.t),self.initialSet)
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
        #print("min: ",min)
        a1=float(nstr(reachSet[ind][0]).split(',')[0][1:])
        b1=float(nstr(reachSet[ind][0]).split(',')[1][:-1])
        a2=self.unsafe[ind][0][0]
        b2=self.unsafe[ind][0][1]
        c=abs(b1-a1)/2
        rad=(b1-c)+1
        #print(rad,c,ind)
        if a2>b1:
            return abs(a2-c)/rad
        elif b2<a1:
            return abs(b2-c)/rad

    def getRobustMetric(self):
        beta=self.getSafeBloat()
        if self.method.lower()=='kagstrom1':
            return self.getRobustMetricKagstrom1(beta)
        elif self.method.lower()=='kagstrom2':
            return self.getRobustMetricKagstrom2(beta)
        elif self.method.lower()=='loan':
            return self.getRobustMetricLoan(beta)

    def getRobustMetricKagstrom1(self,beta):
        kag1=BloatKagstrom(self.A,[])
        (Q,lam,M)=kag1.decompose()
        normM=IntervalNorm.spectralNorm(M)
        Mfac=BloatKagstrom.computeP(normM*self.t,self.A.shape[0])
        robMet=math.log((beta/Mfac)+1)/(Mfac*self.t)
        return robMet

    def getRobustMetricKagstrom2(self,beta):
        kag2=BloatKagstrom(self.A,[])
        (S,N,l,ep)=BloatKagstrom.JNFDecomp(self.A)
        D=BloatKagstrom.getD(l,ep,N)
        KSD=BloatKagstrom.computeK(np.matmul(S,D))
        e_ep=np.exp(ep*self.t)
        robMet=math.log((beta/(KSD*e_ep))+1)/KSD*self.t
        return robMet

    def getRobustMetricLoan(self,beta):
        loan=BloatLoan(self.A,{})
        e_pow=np.exp(math.log(beta)-(loan.computeMu()*self.t)+(loan.computeAlpha()*self.t))
        robMet=SP.lambertw(e_pow.real)/self.t
        return robMet.real

    def getMaxSafeBloat(self):
        beta=self.getSafeBloat()
        reachSet=np.matmul(SLA.expm(self.A*self.t),self.initialSet)
        print(SLA.expm(self.A*self.t))
        print("Reachable Set (Without perturbation)")
        print(reachSet)
        print("-----")
        print("Unsafe Set")
        print(self.unsafe)
        print("-----")
        print("Safe Bloating Factor: ",beta)
        print("Started...")
        b1=self.getRobustMetricKagstrom1(beta)
        b2=self.getRobustMetricKagstrom2(beta)
        b3=self.getRobustMetricLoan(beta)
        print("...Done")
        b=max(b1,b2,b3)
        if (b1==b):
            return (b,'Kagstrom1')
        elif (b2==b):
            return (b,'Kagstrom2')
        else:
            return (b,'Loan')

# Tester
if False:
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

    rM=RobustMetric(A,[],U,method,IS,t,'.',0)
    #print(rM.getRobustMetric())
    print(rM.getMaxSafeBloat())
