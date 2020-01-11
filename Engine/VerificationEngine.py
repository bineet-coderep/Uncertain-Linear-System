'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

from BloatAPI import *
from ReachSetAPI import *
from EigenDecomposition import *
import matplotlib.pyplot as plt

class VerifyBloat:
    '''
    This class provides all the neccessary APIs to compute
    reachable set of a Perturbed Linear Dynamical System
    based on bloating approaches
    '''

    def __init__(self,matA,matB,matE,initset,t,U,m):
        self.A=matA # Matrix A of the Dynamical System
        self.B=matB # Matrix A of the Dynamical System
        self.A_comp=VerifyBloat.createMatrix(self.A,self.B,'.')
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
        self.method=m
        '''
        The method to be used for computation of bloating factor.
        method choices: 'Kagstrom', 'Loan'.
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

    def getBloatingFactor(self,p='slow'):
        if self.method.lower()=='kagstrom1':
            bloat=BloatKagstrom(self.A_comp,self.E)
        elif self.method.lower()=='kagstrom2':
            bloat=BloatLoan(self.A_comp,self.E)
        elif self.method.lower()=='loan':
            bloat=BloatLoan(self.A_comp,self.E)

        bloatFactor=bloat.computeBloatingFactor(self.time,p)
        return bloatFactor

    def computeReachSet(self,p='slow'):
        bloatFactor=self.getBloatingFactor(p)
        reach=ReachSetBloat(self.A_comp,self.initSet,self.time,bloatFactor)
        starReachSet=reach.bloatReachSet(reach.starReachSet())
        return starReachSet

    def computePerturbFreeReachSet(self,p='slow'):
        bloatFactor=self.getBloatingFactor(p)
        reach=ReachSetBloat(self.A_comp,self.initSet,self.time,bloatFactor)
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
        result=ReachSetBloat.isSafe(reachSet,self.Unsafe)
        print("Safety Staus: ",result)

    def isSafePerturbFree(self):
        reachSet=self.computePerturbFreeReachSet()
        result=ReachSetBloat.isSafe(reachSet,self.Unsafe)
        print("Safety Staus: ",result)

    def plotTime(self,start,n,step,p='slow'):
        '''
        Plots the bloating factor with time, upto time tBound
        '''

        if self.method.lower()=='kagstrom1':
            bloat=BloatKagstrom(self.A_comp,self.E)
            (plotX,plotY)=bloat.computeBloatingFactorWithTime(start,n,step,p)
        elif self.method.lower()=='kagstrom2':
            bloat=BloatKagstrom(self.A_comp,self.E)
            (plotX,plotY)=bloat.computeBloatingFactor2WithTime(start,n,step,p)
        elif self.method.lower()=='loan':
            bloat=BloatLoan(self.A_comp,self.E)
            (plotX,plotY)=bloat.computeBloatingFactorWithTime(start,n,step,p)


        plt.autoscale(enable=True, axis='both', tight=False)
        plt.plot(plotX,plotY)
        plt.xlabel("Time")
        plt.ylabel("Bloating Factor")
        plt.show()

    def plotTimeCompare(self,start,n,step,methodList,p='slow'):
        '''
        Plots the bloating factor with time, upto time tBound.
        Comparing all the techniques mentioned in methodList
        in the same plot.
        '''

        plotX=[]
        plotY=[]
        i=1
        for m in methodList:
            if m.lower()=='kagstrom1':
                bloat=BloatKagstrom(self.A_comp,self.E)
                (plotX,plotY)=bloat.computeBloatingFactorWithTime(start,n,step,p)
            elif m.lower()=='kagstrom2':
                bloat=BloatKagstrom(self.A_comp,self.E)
                (plotX,plotY)=bloat.computeBloatingFactor2WithTime(start,n,step,p)
            elif m.lower()=='loan':
                bloat=BloatLoan(self.A_comp,self.E)
                (plotX,plotY)=bloat.computeBloatingFactorWithTime(start,n,step,p)
            plotY=list(filter(lambda x: math.inf!=x,plotY))
            plotX=plotX[:len(plotY)]
            plt.plot(plotX,plotY,label=m)
            i=i+1

        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("Bloating Factor")
        plt.legend(loc='best')
        plt.show()

    def comparePace(self):
        '''
        Plots the bloating factor with time, upto time tBound.
        Comparing the fast and slow techniques mentioned in methodList
        in the same plot.
        '''

        print("==== Bloating Factor as Per Pace ====")
        bloat=BloatLoan(self.A_comp,self.E)
        bFast=bloat.intervalNorm('fast')
        #bSlow=bloat.intervalNorm('slow')
        print("Error Matrix Norm (Fast): ",bFast)
        print("Error Matrix Norm (Slow): ",bSlow)
        print("Relative Difference: ",(bFast-bSlow)/bSlow)
        print("----")

class VerifyDecomp:
    '''
    This class provides all the neccessary APIs to compute
    reachable set of a Perturbed Linear Dynamical System
    based on eigenvalues and eigenvectors decomposition
    '''
    #IMPORTANT: Need Change in ReachSet class

    def __init__(self,matA,matB,B,C,P,initset,t,U):
        self.A=matA # Matrix A of the Dynamical System
        self.B=matB # Matrix A of the Dynamical System
        self.A_comp=VerifyDecomp.createMatrix(self.A,self.B,'.')
        self.b=B
        self.C=C
        self.q=P
        '''
        A_pert=A+Bq^TC
        Please refer to corrollary 2.4 for details
        '''
        self.n=self.A_comp.shape[0]
        self.initialSet=initset # Initial Set.
        '''
        Initial Set is represented using a list of following list
        [C,V,r].
        Where C is the center, V is the matrix of basis vectors with
        first column as 0, r is the radius of the sphere.
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

    def computeReachSet(self):
        reach=ReachSetDecomp(self.A_comp,self.initialSet,self.time,self.b,self.q,self.C)
        return reach.reachSet()

    def computeReachSetPertFree(self):
        reach=ReachSetDecomp(self.A_comp,self.initialSet,self.time,self.b,self.q,self.C)
        return reach.reachSetPertFree()

    def isSafe(self):
        reach=ReachSetDecomp(self.A_comp,self.initialSet,self.T,self.b,self.q,self.c)
        result=reach.isSafe(self.Unsafe)
        print("Safety Staus: ",result)

class VerifyInterval:
    '''
    This class provides all the neccessary APIs to compute
    reachable set of a Perturbed Linear Dynamical System
    based on Interval Arithmetic
    '''

    def __init__(self,matA,matB,matE,initset,t,U):
        self.A=matA # Matrix A of the Dynamical System
        self.B=matB # Matrix A of the Dynamical System
        self.A_comp=VerifyBloat.createMatrix(self.A,self.B,'.')
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

    def computeReachSet(self):
        reach=ReachSetInterval(self.A_comp,self.E,self.initialSet,self.Unsafe,self.time)
        reachSet=reach.reachSet()
        return reachSet

    def computePerturbFreeReachSet(self,p='slow'):
        reach=ReachSetInterval(self.A_comp,self.E,self.initSet,self.time)
        reachSet=reach.reachSetPertFree()
        return reachSet

    def isSafe(self):
        print("Under Construction!!")

    def isSafePerturbFree(self):
        print("Under Construction!!")
