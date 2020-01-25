'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System, Perturbation in the Eigenvectors and Eigenvalues,
compute the reachable set using the following techniques:
  - Perturbation free reach set.
  - Reach set using perturbation, using interval arithmetic.
  - Reach set using perturbation, using reverse eigen decomposition method.

Documentation: Not yet available. (TODO)
'''

import numpy as np
import numpy.linalg as LA
import scipy.linalg as SLA
import math
from time import sleep
import sys
from RevEigenAPI import *

TAYLOR_PRECISION=50000

class ReachSet:
    '''
    Computes the reachable set using the above three mentioned method.
    '''

    def __init__(self,A,B,IS,dVec,dVal,t):
        self.A=ReachSet.createMatrix(A,B,'.')
        self.n=self.A.shape[0]
        self.initialSet=IS # Initial Set, represented using a numpy array
        self.pertEigVal=dVal # Perturbation in Eigen Values
        '''
        Following dictionary data-structure has been used to represent
        a error matrix:
        {
            i: [a,b]
        }
        Indicating i-th eigenvalue can pick any value within the range [a,b]
        '''
        self.pertEigVect=dVec # Perturbation in Eigen Vectors
        '''
        Following dictionary data-structure has been used to represent
        a error matrix:
        {
            (i,j): [a,b]
        }
        Indicating array[i][j] can pick any value within the range [a,b]
        '''
        self.time=t

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

    @staticmethod
    def eigVal2eigMat(mat):
        n=mat.shape[0]
        Lambda=np.zeros((n,n))
        for i in range(n):
            Lambda[i][i]=mat[i]
        return Lambda

    @staticmethod
    def compute_eAt(A):
        n=A.shape[0]
        #matA=np.zeros((n,n),dtype=object)
        matA=np.identity(n,dtype=object)
        matA_i=np.copy(A)

        for i in range(1,TAYLOR_PRECISION):
            sys.stdout.write('\r')
            sys.stdout.write("Progress: "+str((i*100)/TAYLOR_PRECISION)+"%")
            sys.stdout.flush()
            fac=math.factorial(i)
            matA=matA+(matA_i/fac)
            matA_i=np.matmul(matA_i,A)
        print("")

        return matA

    @staticmethod
    def intervalExp(a):
        s=0
        for i in range(TAYLOR_PRECISION):
            s=s+(a**i)/factorial(i)
        return s

    @staticmethod
    def compute_eAt2(A):
        n=A.shape[0]
        matA=np.zeros((n,n),dtype=object)
        for i in range(n):
            matA[i][i]=ReachSet.intervalExp(A[i][i])

        return matA

    def reachSetPertFree(self):
        '''
        Returns the reach set without perturbation
        '''
        (eigVal,eigVec)=LA.eig(self.A*self.time)
        expSig=SLA.expm(ReachSet.eigVal2eigMat(eigVal))
        traj=np.matmul(np.matmul(eigVec,expSig),LA.inv(eigVec))
        #traj=SLA.expm(self.A*self.time)
        rs=np.matmul(traj,self.initialSet)
        #print("m#1: ",np.matmul(SLA.expm(self.A*self.time),self.initialSet))
        return rs

    def reachSetInterval(self):
        '''
        Returns the reach set of the perturbed matrix using interval arithmetic
        '''
        ed=RevEigenDecomp(self.A*self.time,self.pertEigVal,self.pertEigVect)
        Er=ed.getUncertainty()[1]
        #print(Er)
        pertA=(self.A*self.time)+Er
        eA=ReachSet.compute_eAt(pertA)
        rs=np.matmul(eA,self.initialSet)
        return rs

    def reachSetRevED(self):
        '''
        Returns the reach set of the perturbed matrix using
        Reverse Eigen Decomposition
        '''
        (eigVal,eigVec)=LA.eig(self.A*self.time)
        pertEVal=RevEigenDecomp.getPertEigValues(eigVal,self.pertEigVal)
        pertEigVect=RevEigenDecomp.getPertEigVectors(eigVec,self.pertEigVect)
        eSt=ReachSet.compute_eAt2(pertEVal)
        traj=np.matmul(np.matmul(pertEigVect,eSt),IntervalMatrix(pertEigVect).inverse())
        rs=np.matmul(traj,self.initialSet)
        return rs

    def printReport(self):
        print("------------Comparing Various Approaches-----------")
        print("\n-------Reach Set (without Perturbation)-------")
        print(self.reachSetPertFree())
        print("-------\n")
        print("\n-------Reach Set (with Perturbation) using Interval Arithmetic-------")
        print(self.reachSetInterval())
        print("-------\n")
        print("\n-------Reach Set (with Perturbation) using Reverse Eigen Decomposition-------")
        print(self.reachSetRevED())
        print("-------\n")
        print("------------End of Report-----------")


if False:
    print("TEST")
    A=np.array([
    [1,0,1],
    [0,-2,0],
    [1,0,1]
    ])
    B=np.array([])
    IS=np.array([
    [1],
    [1],
    [1],
    ])
    t=5
    dVec={
    (0,0):[0.95,1.05],
    (2,2):[0.95,1.05],
    }
    dVal={
    0:[0.95,1.05]
    }
    (eigVal,eigVec)=LA.eig(A)
    print(eigVal)
    print(eigVec)
    q=ReachSet(A,B,IS,dVec,dVal,t)
    q.printReport()
