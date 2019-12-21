'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

Most the functions are based on the paper:
- 'Bounds and Perturbation Bounds for the Matrix Exponential'
by Bo Kagstrom
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
import math


class BloatKagstrom:
    '''Neccessary APIs required to bloat the Reachable
    Set according to the given uncertainty
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
    def spectralNorm(matA):
        # Computes 2-Norm of matrix matA
        return LA.norm(matA,ord=2)

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

    def centerifyE(self):
        '''
        Break the error matrix to E=[Ac-delta,Ac+delta)
        '''
        Ac=np.zeros((self.n,self.n))
        delta=np.zeros((self.n,self.n))
        for key in self.E:
            Ac[key[0]][key[1]]=(self.E[key][0]+self.E[key][1])/2
            delta[key[0]][key[1]]=self.E[key][1]-Ac[key[0]][key[1]]
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

    def intervalNorm(self):
        '''
        Computes the interval norm of
        E based on Theorem 7 of the
        paper 'Norms of Interval Matrices'
        '''
        norm=-9999
        (Ac,delta)=self.centerifyE()
        for i in range(2**self.n):
            y=BloatKagstrom.generateSignBits(i,self.n,1)
            for j in range(2**self.n):
                z=BloatKagstrom.generateSignBits(j,self.n,0)
                tmp=BloatKagstrom.spectralNorm(Ac+(np.matmul(y,z)*delta))
                #print(Ac+(np.matmul(y,z)**delta))
                if tmp>norm:
                    norm=tmp
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


    def computeBloatingFactor(self,t):
        '''
        Computes the Relative Error Bound
        as per 4.14 (Table 4.1) in the paper
        'Bounds and Perturbation Bounds for the Matrix Exponential'
        '''
        (Q,lam,M)=self.decompose()
        normE=self.intervalNorm()
        normM=BloatKagstrom.spectralNorm(M)

        bloatFactor=BloatKagstrom.computeP(normM*t,self.n)*(np.exp(BloatKagstrom.computeP(normM*t,self.n)*normE*t)-1)

        return bloatFactor


    def computeBloatingFactorWithTime(self,start,n,step):
        '''
        Computes the Relative Error Bound
        as per 4.14 (Table 4.1) in the paper
        'Bounds and Perturbation Bounds for the Matrix Exponential'
        with respect to time.
        '''
        (Q,lam,M)=self.decompose()
        normE=self.intervalNorm()
        normM=BloatKagstrom.spectralNorm(M)

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
