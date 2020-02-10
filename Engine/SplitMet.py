'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.

Documentation: Not yet available. (TODO)
'''

import numpy as np
import numpy.linalg as LA
import mpmath as mp
from gurobipy import *
import sys
import time

BIGM=0.001
EPSILON=1e-3

class Split:
    '''
    Computes the reachable set of a given linear discrete dynamical system A
    form the initial set Theta, upto time T
    '''

    def __init__(self,A,E,theta,T):
        self.A=A # Linear system (without perturbation)
        self.Er=E
        '''
        Following dictionary data-structure has been used to represent
        a error matrix:
        {
            (i,j): [a,b]
        }
        Indicating array[i][j] has a perturbation of -(1-a)% to (b-1)%
        For eg: [0.9,1,1] means a perturbation of -10% and +10%
        '''
        self.Theta=theta # The initial set
        self.T=T #Max number of steps
        self.n=A.shape[0]
        self.Ac=self.computeCenter()
        self.methodName="Split"

    def computeCenter(self):
        '''
        Computes the center point matrix of the interval matrix self.A
        '''
        Ac=np.zeros((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                if (i,j) in self.Er:
                    a=float(self.A[i][j]*self.Er[(i,j)][0])
                    b=float(self.A[i][j]*self.Er[(i,j)][1])
                    c=(a+b)/2
                    Ac[i][j]=c
                else:
                    Ac[i][j]=self.A[i][j]
        return Ac

    def computeU_Interval(self,rs):
        '''
        Computes the effect of uncertainty on the reachable set
        using Interval arithmetic
        '''
        A_tilde=self.computeUncertainMat()
        diff=A_tilde-self.Ac
        U=np.matmul(diff,rs)
        return U

    def computeUncertainMat(self):
        '''
        Computes the interval uncertain matrix
        '''
        A_tilde=np.zeros((self.n,self.n),dtype=object)
        for i in range(self.n):
            for j in range(self.n):
                A_tilde[i][j]=self.A[i][j]
        for key in self.Er:
            a=float(self.Er[key][0]*self.A[key[0]][key[1]])
            b=float(self.Er[key][1]*self.A[key[0]][key[1]])
            #print(a,b)
            #print(mp.mpi(float(a),float(b)))
            #A_tilde[key[0]][key[1]]=mp.mpi(float(self.Er[key][0]*self.A[key[0]][key[1]]),float(self.Er[key][1]*self.A[key[0]][key[1]]))
            #A_tilde[key[0]][key[1]]=mp.mpi(self.Er[key][0]*self.A[key[0]][key[1]],self.Er[key][1]*self.A[key[0]][key[1]])
            A_tilde[key[0]][key[1]]=mp.mpi(min(a,b),max(a,b))

        #print(A_tilde)
        return A_tilde

    def computeU(self,rs):
        '''
        Computes the effect of uncertainty on the reachable set
        using optimization.
        '''

        semiDefFlag=False

        model = Model("qp")
        model.setParam( 'OutputFlag', False )
        #model.params.Presolve=0

        # Create Perturbation Variables
        faultVars=[]
        for key in self.Er:
            name="Pert"+str(key)
            faultVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create Reachable Set Variables
        reachVars=[]
        for i in range(self.n):
            name="IS"+str(i)
            reachVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Add the Perturbation Constraints


        model.optimize() # Updating the model, to use some of the functions

        for var in faultVars:
            kN=var.varName
            i=int(kN.split(',')[0][5:])
            j=int(kN.split(',')[1][:-1])
            k=(i,j)
            name="Pert-C"+str(k)
            model.addConstr(var>=self.Er[k][0],name+".1")
            model.addConstr(var<=self.Er[k][1],name+".2")
        #---------------------------------

        # Add Initial Set Constraints
        for i in range(self.n):
            name="ReachSet-C-"+str(i)
            if (isinstance(rs[i][0],int)) or (isinstance(rs[i][0],np.int64)) or (isinstance(rs[i][0],float)) or (isinstance(rs[i][0],np.float128)):
                model.addConstr(reachVars[i]==rs[i][0],name)
            else:
                a=float(mp.nstr(rs[i][0]).split(',')[0][1:])
                b=float(mp.nstr(rs[i][0]).split(',')[1][:-1])
                if a==b:
                    model.addConstr(reachVars[i]==a,name)
                else:
                    model.addConstr(reachVars[i]>=min(a,b),name+".1")
                    model.addConstr(reachVars[i]<=max(a,b),name+".2")
        #---------------------------------

        # Prepare the objective functions and get min max

        U=np.zeros((self.n,1),dtype=object)

        for i in range(self.n):
            obj=0
            for j in range(self.n):
                if (i,j) in self.Er:
                    pertV=model.getVarByName("Pert"+str((i,j)))
                    obj=obj+((self.A[i][j]*pertV)-self.Ac[i][j])
                else:
                    obj=obj+(self.A[i][j]-self.Ac[i][j])
            obj=obj*reachVars[i]

            # Obtain Minimum
            mn=-9890
            #print(obj," (Min)\n")
            model.setObjective(obj,GRB.MINIMIZE)
            #model.params.Presolve=0
            #model.write("logMin.lp")
            try:
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
                        mn=obj.getValue()
            except:
                semiDefFlag=True


            #print("Min: ",mn)
            #-------------------------------

            # Obtain Maximum
            #print(obj," (Max)\n")
            mx=9890
            model.setObjective(obj,GRB.MAXIMIZE)
            #model.params.Presolve=0
            #model.write("logMax.lp")
            try:
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
                        mx=obj.getValue()
            except:
                #print("Caught it!!")
                semiDefFlag=True


            #print("Max: ",mx)
            #----------------------------------

            #U[i][0]=mp.mpi(mn,mx)
            U[i][0]=(mn,mx)
            #print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ",i)

        if semiDefFlag==True:
            UInterval=self.computeU_Interval(rs)
            for i in range(self.n):
                lb=U[i][0][0]
                ub=U[i][0][1]
                if lb==-9890:
                    lb=np.float(mp.nstr(UInterval[i][0]).split(',')[0][1:])
                if ub==9890:
                    ub=np.float(mp.nstr(UInterval[i][0]).split(',')[1][:-1])
                #print(lb,ub)
                U[i][0]=mp.mpi(lb,ub)
        else:
            for i in range(self.n):
                lb=U[i][0][0]
                ub=U[i][0][1]
                #print(lb,ub)
                U[i][0]=mp.mpi(lb,ub)
        #print("Ret>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        return U

    def getReachableSet(self):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        ORS=self.Theta
        U=self.computeU(ORS)
        t=1
        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Progress: "+str((t*100)/self.T)+"%")
            sys.stdout.flush()
            ORS=np.matmul(self.Ac,ORS)+U
            U=self.computeU(ORS)
            t=t+1
        print()
        return ORS

    def printReachableSet(self):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        start_time=time.time()
        ORS=self.Theta
        U=self.computeU(ORS)
        t=1
        print("-----------------")
        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()
            #print()
            ORS=np.matmul(self.Ac,ORS)+U
            U=self.computeU(ORS)
            t=t+1
        time_taken=time.time()-start_time
        print()

        start_time2=time.time()
        ORSI=self.Theta
        UI=self.computeU_Interval(ORS)
        t=1
        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Interval): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()
            #print()
            ORSI=np.matmul(self.Ac,ORSI)+UI
            UI=self.computeU_Interval(ORSI)
            t=t+1
        time_taken2=time.time()-start_time2
        print()
        print()
        print("\n-------------Reachable Set of the Perturbed System using Splitting Method-------------\n")
        print()
        print("----Optimization----")
        print("Time Taken: ",time_taken)
        print(ORS)
        print("--------")
        print()
        print("----Interval----")
        print("Time Taken: ",time_taken2)
        print(ORSI)
        print("--------")
        print()
        print("Step: ",self.T)
        print("---------------------------------------------------------------")

if False:
    A2=np.array([
    [1,1,-2],
    [2,0.2,0],
    [0,0.1,0.1]
    ])
    A=np.array([
    [2,1,-2,1],
    [2,1.2,0,0],
    [0,0.1,1.1,1],
    [0,0,0,1]
    ])
    E={
    (0,1):[0.9,1.1],
    (1,0):[0.8,1.2],
    (2,3):[0.9,1.1]
    }
    E2={
    (0,1):[0.9,1.5],
    (1,0):[0.8,1.2],
    }
    T=200
    IS2=np.array([
    [1],
    [1],
    [1]
    ])
    IS=np.array([
    [1],
    [1],
    [1],
    [1]
    ])
    Split(A,E,IS,T).printReachableSet()
