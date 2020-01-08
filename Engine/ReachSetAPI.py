'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.
- Class ReachSet: Using Eigenvalue and vector decomposition
- Class ReachSetBloat: The Reachable Set is represented using Generalized Stars according to
the paper 'Parsimonious, Simulation based Verification of Linear Systems'
by Parasara Sridhar Duggirala and Mahesh Viswanathan.

Documentation: Not yet available. (TODO)
'''


import numpy as np
import scipy.linalg as LA
from gurobipy import *
from EigenDecomposition import *

TAYLOR_PRECISION=100


class ReachSetBloat:
    '''
    Compute Reachable Set and bloat it with a given factor

    NOTE: The 'shape' of the reach set is hardcoded to Sphere
    '''
    def __init__(self, matA, initSet, t, blt):
        self.matA=matA # Matrix A of size n*n, represented as a numpy array.
        self.n=self.matA.shape[0] # Dimension of the System
        self.initialSet=initSet # Initial Set.
        '''
        Initial Set is represented using a list of following list
        [C,V,r].
        Where C is the center, V is the matrix of basis vectors with
        first column as 0, r is the radius of the sphere.
        '''
        self.T=t # Time T, upto which ReachSet is to be calculated
        self.bloatingFactor=blt # Bloat the Reachable Set to factor bloatingFactor
        self.radius=r # The radius of the initial set, represented as a sphere

    def compute_eAt(self,vi):
        #vi indicates the vector
        eAt=np.matmul((LA.expm(self.matA*self.T)),self.initialSet[0].reshape(-1,1)+self.initialSet[1][:,vi].reshape(-1,1)) # Computes e^{At}
        return eAt

    def starReachSet(self):
        '''
        Computes the Reachable Set based on Star Representation
        '''
        V_reach=np.zeros(self.initialSet[1].shape) # Basis vectors of the Reachable Set
        r_reach=self.initialSet[2] # Radius of the sphere representing the Reachable Set
        C_reach=self.compute_eAt(0) # Center of the new Reachable Set

        # Algorithm 1 of the above mentioned paper
        for i in range(1,self.n+1):
            tmp=self.compute_eAt(i)
            t=tmp-C_reach
            for j in range(self.n):
                V_reach[j][i]=t[j][0]

        return [C_reach,V_reach,r_reach]

    def bloatReachSet(self,starSet):
        starSet[3]=starSet[3]+self.bloatingFactor
        return starSet

    @staticmethod
    def isSafe(starSet,U):

        '''
        Unsafe Condition:
            startSet>=U
        '''

        model = Model("qp")

        C=starSet[0] # Center
        V=starSet[1][:,1:] # Basis Vectors
        r=starSet[2] # Radius of the sphere

        alpha=[]
        # Define the variables alpha
        for i in range(self.n):
            alpha.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name="alpha"+str(i),vtype='C'))
        #---------------------------

        # Add constraints related to alpha
        alphaC=[]
        for i in range(self.n):
            constr=C[i]
            for j in range(self.n):
                constr=constr+(alpha[i]*V[i][j])
            alphaC.append(constr)

        for i in range(self.n):
            nam="AlphaC-"+str(i)+"-"
            model.addConstr(alphaC[i]>=U[i])
        #---------------------------------

        # Add Predicate Conditions
        for i in range(self.n):
            nam="PredicateC-"+str(i)+"-"
            constr=((alpha[i]-C[i])*(alpha[i]-C[i]))
            model.addConstr(constr<=r)
        #----------------------------------


        # Solve the model
        model.setObjective(1)
        model.optimize()
        status = model.Status
        if status==GRB.Status.UNBOUNDED:
            print("UNBOUNDED ")
        else:
            if status == GRB.Status.INF_OR_UNBD or \
               status == GRB.Status.INFEASIBLE  or \
               status == GRB.Status.UNBOUNDED:
                print('**The model cannot be solved because it is infeasible or unbounded**')
            else:
                print("\n\nCounter Example\n\n")
                for v in model.getVars():
                    print('%s %g' % (v.varName, v.x))
                return False


            print("------------------------------")
        #------------------------------------------------------------------------------------

        return True

class ReachSetDecomp:
    '''
    Compute Reachable Set using Eigenvalue and Eigenvector
    decomposition

    '''
    def __init__(self, A, initSet, t, B, P, C):
        '''
        A_pert=A+BPC
        Please refer to corrollary 2.4 for details
        '''
        self.A=A
        self.b=B
        self.C=C
        self.q=P
        self.n=A.shape[0] # Dimension of the System
        self.initialSet=initSet # Initial Set.
        '''
        Initial Set is represented using a list of following list
        as a row vector (numpy array)
        '''
        self.T=t # Time T, upto which ReachSet is to be calculated

    @staticmethod
    def computeSigEAt(matA,t):
        l=matA.shape[0]
        A=np.zeros((l,l),dtype=object)
        for i in range(l):
            A[i][i]=matA[i]
        A=A*t
        S=np.zeros((l,l),dtype=object)
        for i in range(l):
            S[i][i]=ReachSetDecomp.intervalExp(A[i][i])
        #print("S: ",S)
        return S

    @staticmethod
    def intervalExp(a):
        s=0
        for i in range(TAYLOR_PRECISION):
            s=s+(a**i)/factorial(i)
        return s

    def reachSet(self):
        (eval,evects)=EigenDecompose(self.A,self.b,self.q,self.C).decompose()
        print("eval: ",eval)
        eSig=ReachSetDecomp.computeSigEAt(eval,self.T)
        #print("eSig: ",eSig)
        evectsInv=IntervalMatrix(evects).inverse()
        rSet=np.matmul(np.matmul(evects,eSig),evectsInv)
        rSet=np.matmul(rSet,self.initialSet)
        return rSet

    def reachSetPertFree(self):
        (eval,evects)=LA.eig(self.A)
        print("eval: ",eval)
        print("evects: \n",evects)
        l=self.A.shape[0]
        S=np.zeros((l,l),dtype=object)
        for i in range(l):
            S[i][i]=np.exp(eval[i])
        #print("S: ",S )
        rSet=np.matmul(np.matmul(evects,S),LA.inv(evects))
        rSet=np.matmul(rSet,self.initialSet)
        return rSet

    def isSafe(self,U):
        range=np.matmul(self.reachSet,self.initialSet)
        '''
        Unsafe Condition:
            startSet>=U
        '''
        for i in range(self.n):
            if float(nstr(range[i]).split(',')[1][:-1]) >= U[i]:
                return False

        return True
