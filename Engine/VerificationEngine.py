'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

from BloatAPI import *
from ReachSetAPI import *
import matplotlib.pyplot as plt

class Verify:
    '''
    This class provides all the neccessary APIs to compute
    reachable set of a Perturbed Linear Dynamical System
    '''

    def __init__(self,matA,matB,matE,initset,t,U):
        self.A=matA # Matrix A of the Dynamical System
        self.B=matB # Matrix A of the Dynamical System
        self.A_comp=Verify.createMatrix(self.A,self.B,'.')
        self.E=matE # Matrix E, represents error.
        '''
        Following dictionary data-structure has been used to represent
        the error matrix E:
        {
            (i,j): [a,b]
        }
        Indicating E[i][j] can pick any value within the range [a,b]
        '''
        self.n=self.A_comp.shape[0]
        self.initialSet=initset # Initial Set.
        '''
        Initial Set is represented using a list of following list
        [C,V,r].
        Where C is the center, V is the matrix of basis vectors with
        first column as 0, r is the radius of the sphere.
        '''
        self.E=matE # Matrix E, represents error.
        '''
        Following dictionary data-structure has been used to represent
        the error matrix E:
        {
            (i,j): [a,b]
        }
        Indicating E[i][j] can pick any value within the range [a,b]
        '''
        self.n=0 # Dimension of the System
        self.time=t # Time upto which Verification is to be performed
        self.Unsafe=U # Unsafe set, list of numbers
        '''
        The unsafe condition for a state i is the following:
        state of i >= U[i]
        '''

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

        C=np.zeros((n,n),dtype=np.float128)
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

    def getBloatingFactor(self):
        bloat=Bloat(self.A_comp,self.E)
        bloatFactor=bloat.computeBloatingFactor(self.time)
        return bloatFactor

    def computeReachSet(self):
        bloatFactor=self.getBloatingFactor()
        reach=ReachSet(self.A_comp,self.initSet,self.time,bloatFactor)
        starReachSet=reach.bloatReachSet(reach.starReachSet())
        return starReachSet

    def computePerturbFreeReachSet(self):
        bloatFactor=self.getBloatingFactor()
        reach=ReachSet(self.A_comp,self.initSet,self.time,bloatFactor)
        starReachSet=reach.starReachSet()
        return starReachSet

    def printReachSet(self,reachSet):
        print("======== Reach Set =========")
        print("\t ----Center----")
        print(reachSet[0])
        print("\t ----Basis----")
        print(reachSet[1])
        print("\t ----Radius----")
        print(reachSet[2])

    def isSafe(self):
        reachSet=self.computeReachSet()
        result=ReachSet.isSafe(reachSet,self.Unsafe)
        print("Safety Staus: ",result)

    def isSafePerturbFree(self):
        reachSet=self.computePerturbFreeReachSet()
        result=ReachSet.isSafe(reachSet,self.Unsafe)
        print("Safety Staus: ",result)

    def plotTime(self,start,n,step):
        '''
        Plots the bloating factor with time, upto time tBound
        '''

        bloat=Bloat(self.A_comp,self.E)
        #print(bloat.computeBloatingFactor(0.001))
        (plotX,plotY)=bloat.computeBloatingFactorWithTime(start,n,step)

        plt.autoscale(enable=True, axis='both', tight=False)
        plt.plot(plotX,plotY)
        plt.xlabel("Time")
        plt.ylabel("Bloating Factor")
        plt.show()
