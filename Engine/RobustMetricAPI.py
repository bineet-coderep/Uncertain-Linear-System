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

class RobustMetric:
    '''
    Computes Robustness Metric according to the following
    three Bloating Methods:
        1. Kagstrom 1
        2. Kagstrom 2
        3. Loan
    '''

    def __init__(self,A,U,m,t):
        self.A=A # The Linear Dynamical Matrix A
        self.unsafe=U # The Unsafe Set
        self.method=m # Method to be used to compute Robustness Metric (Kagstrom1, Kagstrom2, Loan)
        self.t=t # Time

    def getSafeBloat(self):
        print("Under Construction!!")
        beta=1000
        return beta

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


# Tester
if True:
    A=np.array([
    [0,0,1,0,1],
    [1,0,0,1,1],
    [0,0,0,0,0],
    [0,0,1,0,1],
    [0,1,1,0,0],
    ])
    U=[]
    t=20
    method='Kagstrom2'

    rM=RobustMetric(A,U,method,t)
    print(rM.getRobustMetric())
