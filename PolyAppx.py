'''
Author: Bineet Ghosh
Supervised By: Prof. Sridhar
Email: ghosh.bineet22@gmail.com
'''

'''
This code tests an idea where the polynomials
representing the uncertainties are approximated.
'''



import numpy as np
import itertools
import time
import math
from gurobipy import *
from numpy import linalg as LA
import xml.etree.ElementTree as ET

M=1e3

class Polynomial:
    '''This class represents a reference Polynomial
    and approximates polynomials based on the
    reference polynomial
    '''

    '''
    Ref Polynomial: x^4+y^4+z^4+xyz
    (Hardcoded)
    '''

    def __init__(self,poly,b,si):
        self.actualPoly=poly
        self.bounds=b
        self.appxSize=4
        self.size=10
        self.sampleInterval=si
        self.appxPoly=self.approximatePolyCheby()

    def polyMul(self,P):
        pl=[]
        for i in range(self.appxSize):
            tmp=[]
            for j in range(self.appxSize):
                tmp.append(self.appxPoly[i]*P.appxPoly[j])
            pl=pl+tmp

        #Hardcoded
        pl[1]=pl[1]+pl[4]
        pl[2]=pl[2]+pl[8]
        pl[3]=pl[3]+pl[12]
        pl[4]=pl[5]
        pl[5]=pl[6]+pl[9]
        pl[6]=pl[7]+pl[13]
        pl[7]=pl[11]+pl[14]
        pl[8]=pl[10]
        pl[9]=pl[15]

        return Polynomial(pl,self.bounds,self.sampleInterval)

    def polyPower(self,p):
        pt=Polynomial(self.actualPoly,self.bounds,self.sampleInterval)
        for i in range(p-1):
            pt=self.polyMul(pt)
        return pt

    def approximatePolyOpti(self):

        #Hardcoded

        flag=True
        model = Model("qp")

        #Adding polynomial variables
        x=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="x",vtype='C')
        y=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="y",vtype='C')
        z=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="z",vtype='C')

        #Reference Polynomial Variables
        refPolyVars=[]
        for i in range(self.appxSize):
            name="Alpha ("+str(i)+")"
            refPolyVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))

        #Epsilon Variable
        epsilon=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="epsilon",vtype='C')

        #Adding Bounds to Polynomial Variables
        #model.addConstr(x==self.bounds[0][0],"x.lb")
        model.addConstr(x==self.bounds[0][1],"x.ub")
        #model.addConstr(y==self.bounds[1][0],"y.lb")
        model.addConstr(y==self.bounds[1][1],"y.ub")
        #model.addConstr(z==self.bounds[2][0],"z.lb")
        model.addConstr(z==self.bounds[2][1],"z.ub")

        #Adding Approximation Constraints
        model.addConstr(refPolyVars[0]-self.actualPoly[0]<=epsilon)
        model.addConstr(refPolyVars[0]-self.actualPoly[0]>=-epsilon)
        model.addConstr(refPolyVars[1]-self.actualPoly[4]<=epsilon)
        model.addConstr(refPolyVars[1]-self.actualPoly[4]>=-epsilon)
        model.addConstr(refPolyVars[2]-self.actualPoly[8]<=epsilon)
        model.addConstr(refPolyVars[2]-self.actualPoly[8]>=-epsilon)
        model.addConstr(refPolyVars[3]-(self.actualPoly[1]+self.actualPoly[2]+self.actualPoly[3]+self.actualPoly[5]+self.actualPoly[6]+self.actualPoly[7]+self.actualPoly[9])<=epsilon)
        model.addConstr(refPolyVars[3]-(self.actualPoly[1]+self.actualPoly[2]+self.actualPoly[3]+self.actualPoly[5]+self.actualPoly[6]+self.actualPoly[7]+self.actualPoly[9])>=-epsilon)

        #Adding Bounds to Ref Polynomial
        for i in range(self.appxSize):
            model.addConstr(refPolyVars[i]<=M)
            model.addConstr(refPolyVars[i]>=-M)

        #Objective Function
        obj=epsilon

        model.setObjective(obj,GRB.MINIMIZE)
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
                #print('Obj: %g' % obj.getValue())
                return True


            print("------------------------------")

        return False

    def approximatePolyCheby(self):

        points=self.sampler(self.sampleInterval)

        A=np.zeros((len(points),self.appxSize))
        B=np.zeros((len(points),1))

        for i in range(len(points)):
            x=points[i][0]
            y=points[i][1]
            z=points[i][2]
            A[i][0]=x*x*x*x
            A[i][1]=y*y*y*y
            A[i][2]=z*z*z*z
            A[i][3]=x*y*z

            B[i][0]=(np.power(x,8)*self.actualPoly[0])+((np.power(x,4)*np.power(y,4))*self.actualPoly[1])+((np.power(x,4)*np.power(z,4))*self.actualPoly[2])+((np.power(x,5)*y*z)*self.actualPoly[3])+((np.power(y,8))*self.actualPoly[4])+((np.power(y,4)*np.power(z,4))*self.actualPoly[5])+((x*np.power(y,5)*z)*self.actualPoly[6])+((x*y*np.power(z,5))*self.actualPoly[7])+((np.power(z,8))*self.actualPoly[8])+((np.power(x,2)*np.power(y,2)*np.power(z,2))*self.actualPoly[9])
        AI=LA.pinv(A)
        X=np.matmul(AI,B)

        return (list(X.reshape(1,self.appxSize)[0]))

    def sampler(self,g):
        stepX=(self.bounds[0][1]-self.bounds[0][0])/g
        stepY=(self.bounds[1][1]-self.bounds[1][0])/g
        stepZ=(self.bounds[2][1]-self.bounds[2][0])/g

        initX=self.bounds[0][0]
        initY=self.bounds[1][0]
        initZ=self.bounds[2][0]

        pts=[]
        for i in range(g*g*g):
            pts.append((initX,initY,initZ))
            initX=initX+stepX
            initY=initY+stepY
            initZ=initZ+stepZ

        return pts

    def printPoly(self):
        print(self.actualPoly[0]," x^4 + ",self.actualPoly[1]," x^4.y^4 + ",self.actualPoly[2]," x^4.z^4 + ",self.actualPoly[3]," x^5.y.z + ",self.actualPoly[4]," y^8 + ",self.actualPoly[5]," y^4.z^4 + ",self.actualPoly[6]," x.y^5.z + ",self.actualPoly[7]," x.y.z^5 + ",self.actualPoly[8]," z^8 + ", self.actualPoly[9],"+ x^2.y^2.z^2")

    def printAppxPoly(self):
        print(self.appxPoly[0]," x^4 + ",self.appxPoly[1]," y^4 + ",self.appxPoly[2]," z^4 + ",self.appxPoly[3]," x.y.z ")


class Verification:
    '''This is to test the scalibility of the method
    '''

    def __init__(self,poly,t):
        self.polynomial=poly
        self.steps=t

    def isSafe(self):

        P=self.polynomial
        for i in range(self.steps):
            print("=========== Time ",i,"============")
            if (self.isSafeTime(P)!=True):
                return False
        return True


    def isSafeTime(self,P):
        flag=False
        model = Model("qp")

        #Adding polynomial variables
        x=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="x",vtype='C')
        y=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="y",vtype='C')
        z=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="z",vtype='C')

        #Adding Bounds to Polynomial Variables
        model.addConstr(x>=P.bounds[0][0],"x.lb")
        model.addConstr(x<=P.bounds[0][1],"x.ub")
        model.addConstr(y>=P.bounds[1][0],"y.lb")
        model.addConstr(y<=P.bounds[1][1],"y.ub")
        model.addConstr(z>=P.bounds[2][0],"z.lb")
        model.addConstr(z<=P.bounds[2][1],"z.ub")

        #Encoding the Approximated Polynomial
        polyAppx=(P.appxPoly[0]*(x*x))+(P.appxPoly[1]*(y*y))+(P.appxPoly[2]*(z*z))+(P.appxPoly[3]*(x*y))
        model.addConstr(polyAppx>=7)

        model.setObjective(x)
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
                #print('Obj: %g' % obj.getValue())
                return False

            P=P.polyMul(self.polynomial)
            print("------------------------------")

            return True






p1=Polynomial([1,0,1,0,1,0,0,0.8,1,0.2],[[-1,1],[-1,1],[-1,1]],4)
p2=Polynomial([0,1,0,1,0,1,1,0,1,0.6],[[-1,1],[0,1],[-1,0]],4)

v1=Verification(p1,20)
print(v1.isSafe())
