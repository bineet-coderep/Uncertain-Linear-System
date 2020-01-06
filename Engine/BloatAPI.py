'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

Most the functions are based on the paper:
- 'Bounds and Perturbation Bounds for the Matrix Exponential'
by Bo Kagstrom
- 'The Sensitivity of the Matrix Exponential'
by Chales Van Loan
- 'Norms of Interval Matrices'
by Raena Farhadsefat, Jirı Rohn and Taher Lotf


- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.
- Given a perturbation, the bloating factor according to which e^{At} is to be
bloated to accomodate e^{(A+B)t} is computed.

Documentation: Not yet available. (TODO)
'''


import numpy as np
import numpy.linalg as LA
import scipy.linalg as SLA
from sympy.matrices import *
import math
from random import seed
from random import random

class BloatKagstrom:
    '''
    Neccessary APIs required to bloat the Reachable
    Set according to the given uncertainty. Based on the paper
    'Bounds and Perturbation Bounds for the Matrix Exponential'
    by Bo Kagstrom.
    '''
    def __init__(self, matA, matE):
        self.A=matA # Matrix A of size n*n, represented as a numpy array.
        self.E=matE # Matrix E, represents error.
        '''
        Following dictionary data-structure has been used to represent
        the error matrix E:
        {
            (i,j): [a,b]
        }
        Indicating E[i][j] can pick any value within the range [a,b]
        '''
        self.n=self.A.shape[0] # Dimension of the System

    @staticmethod
    def computeP(x,n):
        '''
        compute p_{n-1}(x) according to the paper
        'Bounds and Perturbation Bounds for the Matrix Exponential'
        '''
        s=0
        for k in range(n-1):
            xk=x**k
            k_fact=math.factorial(k)
            s=s+(xk/k_fact)
        return s

    def intervalNorm(self,p='slow'):
        '''
        Computes the interval norm of
        E based on Theorem 7/10 of the
        paper 'Norms of Interval Matrices'
        '''
        norm=IntervalNorm(self.E,self.n,p).getNorm()
        return norm

    def decompose(self):
        '''
        Decompose A=QTQ^H
        Q: Unitary
        T=lam+M
        lam: Diagonal
        M: Strict Upper Triangular

        This functions returns: (Q,lam,M)

        IMPORTANT: The M is not minimal
        norm in this implementation
        '''

        (T,Q)=SLA.schur(self.A)
        lam=np.zeros((self.n,self.n))
        for i in range(self.n):
            lam[i][i]=T[i][i]
            T[i][i]=0

        return (Q,lam,T)

    def computeBloatingFactor(self,t,p='slow'):
        '''
        Computes the Relative Error Bound
        as per 4.14 (Table 4.1) in the paper
        'Bounds and Perturbation Bounds for the Matrix Exponential'
        '''
        (Q,lam,M)=self.decompose()
        normE=self.intervalNorm(p)
        normM=IntervalNorm.spectralNorm(M)

        bloatFactor=BloatKagstrom.computeP(normM*t,self.n)*(np.exp(BloatKagstrom.computeP(normM*t,self.n)*normE*t)-1)

        return bloatFactor

    def computeBloatingFactorWithTime(self,start,n,step,p='slow'):
        '''
        Computes the Relative Error Bound
        as per 4.14 (Table 4.1) in the paper
        'Bounds and Perturbation Bounds for the Matrix Exponential'
        with respect to time.
        '''

        (Q,lam,M)=self.decompose()
        normE=self.intervalNorm(p)
        normM=IntervalNorm.spectralNorm(M)

        print("Norm of E: ",normE)
        print("Norm of M: ", normM)
        print("")

        timeAxis=[]
        fAxis=[]
        t=start
        it=0
        while True:
            bloatFactor=BloatKagstrom.computeP(normM*t,self.n)*(np.exp(BloatKagstrom.computeP(normM*t,self.n)*normE*t)-1)
            timeAxis.append(t)
            print("Time ",t,": ",bloatFactor)
            fAxis.append(bloatFactor)
            t=t+step
            it=it+1
            if (it>n):
                break

        return (timeAxis,fAxis)

    def computeBloatingFactor2(self,t,p='slow'):
        '''
        Computes the Relative Error Bound
        as per 4.12 (Table 4.1) in the paper
        'Bounds and Perturbation Bounds for the Matrix Exponential'
        '''

        (S,N,l,ep)=BloatKagstrom.JNFDecomp(self.A)
        D=BloatKagstrom.getD(l,ep,N)
        K=BloatKagstrom.computeK(np.matmul(S,D))
        normE=self.intervalNorm(p)
        bloat=K*np.exp(ep*t)*(np.exp(K*normE*t)-1)
        return bloat


    def computeBloatingFactor2WithTime(self,start,n,step,p='slow'):
        '''
        Computes the Relative Error Bound
        as per 4.12 (Table 4.1) in the paper
        'Bounds and Perturbation Bounds for the Matrix Exponential'
        with respect to time.
        '''


        (S,N,l,ep)=BloatKagstrom.JNFDecomp(self.A)
        D=BloatKagstrom.getD(l,ep,N)
        K=BloatKagstrom.computeK(np.matmul(S,D))
        normE=self.intervalNorm(p)

        timeAxis=[]
        fAxis=[]
        t=start
        it=0
        while True:
            bloatFactor=K*np.exp(ep*t)*(np.exp(K*normE*t)-1)
            timeAxis.append(t)
            print("Time ",t,": ",bloatFactor)
            fAxis.append(bloatFactor)
            t=t+step
            it=it+1
            if (it>n):
                break

        return (timeAxis,fAxis)

    @staticmethod
    def computeK(S):
        '''
        Computes K(S) as given in the paper
        The Sensitivity of the Matrix Exponential
        by Kagstrom

        K(S)=||S||*||S^-1||
        '''

        K=IntervalNorm.spectralNorm(S)*(IntervalNorm.spectralNorm(LA.inv(S)))
        return K

    @staticmethod
    def JNFDecomp(A):
        '''
        Decomposes the matrix A to Jordan Normal Form
        A=SJS^-1
        returns S,N,len,ep
        Where N: has 1 in the same positions (i,i+1) as J
        len: countJordanBlocks(J), the length of Jordan Blocks
        ep: 0 < ep < -alpha(A). alpha(A) is the maximum eigen value (real part) [Yet to implement]
        '''

        a=Matrix(A)
        (s,j)=a.jordan_form()
        print("Starting....")
        S=np.array(s)
        J=np.array(j)
        S=S.astype('float')
        J=J.astype('float')
        N=J
        for i in range(N.shape[0]):
            N[i][i]=0

        '''seed(1)
        ep=(100*random())%(-max(LA.eig(A)[0].real))
        if ep<=0:
            print("Something went wrong!!")
            print("Max Eigen value is: ",max(LA.eig(A)[0].real))
            exit(0)
        '''
        ep=0.05
        print("....Done")

        return (S,N,BloatKagstrom.countJordanBlocks(J),ep)

    @staticmethod
    def countJordanBlocks(J):
        '''
        Counts the lengths of each Jordan Blocks
        '''
        n=J.shape[0]
        l=[]
        c=0
        for i in range(n-1):
            if J[i][i+1]!=1:
                l.append(c+1)
                c=0
            else:
                c=c+1
        if sum(l)!=n:
            l.append(1)

        return l

    @staticmethod
    def getD(leng,epsilon,N):
        '''
        Returns the matrix according to the given length leng
        as given in the paper The Sensitivity of the Matrix Exponential
        by Kagstrom
        '''

        delta=min(1,epsilon/LA.norm(N,ord='fro'))
        n=sum(leng)
        D=np.zeros((n,n))
        ind=0
        for c in leng:
            for j in range(c):
                D[ind][ind]=delta**j
                ind=ind+1

        return D


class BloatLoan:
    '''
    Neccessary APIs required to bloat the Reachable
    Set according to the given uncertainty. Based on the paper
    'The Sensitivity of the Matrix Exponential'
    by Chales Van Loan
    '''

    def __init__(self, matA, matE):
        self.A=matA # Matrix A of size n*n, represented as a numpy array.
        self.E=matE # Matrix E, represents error.
        '''
        Following dictionary data-structure has been used to represent
        the error matrix E:
        {
            (i,j): [a,b]
        }
        Indicating E[i][j] can pick any value within the range [a,b]
        '''
        self.n=self.A.shape[0] # Dimension of the System

    def intervalNorm(self,p='slow'):
        '''
        Computes the interval norm of
        E based on Theorem 7 of the
        paper 'Norms of Interval Matrices'
        '''
        norm=IntervalNorm(self.E,self.n,p).getNorm()
        return norm

    def computeBloatingFactor(self,t,p='slow'):
        '''
        Computes the Relative Error Bound
        as per Theorem 2 in the paper
        'The Sensitivity of the Matrix Exponential'
        '''
        normE=self.intervalNorm(p)
        alphaA=self.computeAlpha()
        muA=self.computeMu()

        ePow=(muA-alphaA+normE)*t
        bloatFactor=t*normE*np.exp(ePow)

        return bloatFactor

    def computeBloatingFactorWithTime(self,start,n,step,p='slow'):
        '''
        Computes the Relative Error Bound
        as per Theorem 2 in the paper
        'The Sensitivity of the Matrix Exponential'
        with respect to time.
        '''
        normE=self.intervalNorm(p)
        alphaA=self.computeAlpha()
        muA=self.computeMu()
        ePow=(muA-alphaA+normE)

        timeAxis=[]
        fAxis=[]
        t=start
        it=0
        while True:
            bloatFactor=t*normE*np.exp(ePow*t)
            timeAxis.append(t)
            print("Time ",t,": ",bloatFactor)
            fAxis.append(bloatFactor)
            t=t+step
            it=it+1
            if (it>n):
                break

        return (timeAxis,fAxis)

    def computeAlpha(self):
        '''
        Computes alpha(A) as per the paper 'The Sensitivity of the Matrix Exponential'
        by Chales Van Loan
        '''
        return max(LA.eig(self.A)[0].real)

    def conjugateFactor(self):
        '''
        Returns the following:
        (A*+(A/2))
        '''
        AStar=self.A.conjugate().transpose()
        return (AStar+(self.A/2))

    def computeMu(self):
        '''
        Computes mu(A) as per the paper 'The Sensitivity of the Matrix Exponential'
        by Chales Van Loan
        '''
        return (max(LA.eig(self.conjugateFactor())[0]))


class IntervalNorm:
    '''
    Computes Norm-2 of a Interval Matrix
    based on the paper 'Norms of Interval Matrices'
    by Raena Farhadsefat, Jirı Rohn and Taher Lotf
    '''

    def __init__(self,matrix,s,p='slow'):
        self.A=matrix
        self.pace=p
        self.n=s

    def getNorm(self):
        '''
        Return: Norm2 by default,
        Frobenius Norm if fast method is wanted
        '''
        if self.pace.lower() == 'slow':
            return self.intervalNorm2()
        else:
            return self.frobeniusNorm()

    def centerify(self):
        '''
        Break the error matrix to A=[Ac-delta,Ac+delta)
        '''
        Ac=np.zeros((self.n,self.n))
        delta=np.zeros((self.n,self.n))
        for key in self.A:
            Ac[key[0]][key[1]]=(self.A[key][0]+self.A[key][1])/2
            delta[key[0]][key[1]]=self.A[key][1]-Ac[key[0]][key[1]]
        return (Ac,delta)

    @staticmethod
    def generateSignBits(n,size,axis):
        '''
        generates a list of (+,-1) of size n,
        based on n's binary interpretation
        '''
        s=np.binary_repr(n,size)

        if axis==0:
            bit=np.zeros((1,size))
            for i in range(size):
                if s[i]=='1':
                    bit[0][i]=1
                else:
                    bit[0][i]=-1
            return bit
        else:
            bit=np.zeros((size,1))
            for i in range(size):
                if s[i]=='1':
                    bit[i][0]=1
                else:
                    bit[i][0]=-1
            return bit

    def intervalNorm2(self):
        '''
        Computes the interval norm of
        A based on Theorem 7 of the
        paper 'Norms of Interval Matrices'
        '''
        #print("SLOW")
        norm=-9999
        (Ac,delta)=self.centerify()
        for i in range(2**self.n):
            y=IntervalNorm.generateSignBits(i,self.n,1)
            for j in range(2**self.n):
                z=IntervalNorm.generateSignBits(j,self.n,0)
                tmp=IntervalNorm.spectralNorm(Ac+(np.matmul(y,z)*delta))
                if tmp>norm:
                    norm=tmp
        return norm

    @staticmethod
    def spectralNorm(matA):
        # Computes 2-Norm of matrix matA
        return LA.norm(matA,ord=2)

    def frobeniusNorm(self):
        '''
        Computes the interval norm of
        A based on Theorem 10 of the
        paper 'Norms of Interval Matrices'
        '''
        #print("FAST")
        (Ac,delta)=self.centerify()
        Ac=abs(Ac)
        return LA.norm(Ac+delta,ord='fro')



# Tester --------------

if False:
    A=np.array([
    [1,0,0,1,0],
    [0,1,0,0,1],
    [0,0,1,0,0],
    [1,0,0,1,1],
    [0,0,0,0,1],
    ])
    E={
    (0,2): [-0.2,0.2],
    (3,2): [-0.1,0.1]
    }
    b=BloatKagstrom(A,E)
    print(b.computeBloatingFactor2WithTime(0,5,0.01))
