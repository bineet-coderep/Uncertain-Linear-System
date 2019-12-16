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
    Ref Polynomial: x^2+y^2+z^2+xy+yz+xz
    (Hardcoded)
    '''

    def __init__(self,poly,b,si):
        self.actualPoly=poly
        self.bounds=b
        self.appxSize=6
        self.size=15
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
        pl[1]=pl[1]+pl[6]+pl[21]
        pl[2]=pl[2]+pl[12]+pl[35]
        pl[3]=pl[3]+pl[18]
        pl[4]=pl[4]+pl[23]+pl[33]+pl[24]
        pl[5]=pl[5]+pl[30]
        pl[6]=pl[7]
        pl[7]=pl[8]+pl[13]+pl[28]
        pl[8]=pl[9]+pl[19]
        pl[9]=pl[10]+pl[25]
        pl[10]=pl[11]+pl[22]+pl[27]+pl[31]
        pl[11]=pl[14]
        pl[12]=pl[15]+pl[20]+pl[29]+pl[34]
        pl[13]=pl[16]+pl[26]
        pl[14]=pl[17]+pl[32]

        return Polynomial(pl,self.bounds,self.sampleInterval)

    def polyPower(self,p):
        pt=Polynomial(self.actualPoly,self.bounds,self.sampleInterval)
        for i in range(p-1):
            pt=self.polyMul(pt)
        return pt

    def approximatePolyCheby(self):

        points=self.sampler(self.sampleInterval)

        A=np.zeros((len(points),self.appxSize))
        B=np.zeros((len(points),1))

        for i in range(len(points)):
            x=points[i][0]
            y=points[i][1]
            z=points[i][2]
            A[i][0]=x*x
            A[i][1]=y*y
            A[i][2]=z*z
            A[i][3]=x*y
            A[i][4]=y*z
            A[i][5]=x*z
            sum1=(self.actualPoly[0]*(x**4))+(self.actualPoly[1]*((x**2)*(y**2)))+(self.actualPoly[2]*((x**2)*(z**2)))+(self.actualPoly[3]*((x**3)*(y)))+(self.actualPoly[4]*((x**2)*(y)*(z)))
            sum2=(self.actualPoly[5]*((x**3)*(z)))+(self.actualPoly[6]*((y**4)))+(self.actualPoly[7]*((y**2)*(z**2)))+(self.actualPoly[8]*((x)*(y**3)))+(self.actualPoly[9]*((y**3)*(z)))
            sum3=(self.actualPoly[10]*((x)*(y**2)*(z)))+(self.actualPoly[11]*((z**4)))+(self.actualPoly[12]*((x)*(y)*(z**2)))+(self.actualPoly[13]*((y)*(z**3)))+(self.actualPoly[14]*((x)*(z**3)))
            B[i][0]=sum1+sum2+sum3
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
        print(self.actualPoly[0]," x^4 + ",self.actualPoly[1]," x^2.y^2 + ",self.actualPoly[2]," x^2.z^2 + ",self.actualPoly[3]," x^3.y + ",self.actualPoly[4]," x^2.y.z + ",self.actualPoly[5]," x^3.z + ",self.actualPoly[6]," y^4 + ",self.actualPoly[7]," y^2.z^2 + ",self.actualPoly[8]," x.y^3 + ")
        print(self.actualPoly[9]," y^3.z + ",self.actualPoly[10]," x.y^2.z + ",self.actualPoly[11]," z^4 + ",self.actualPoly[12]," x.y.z^2 + ",self.actualPoly[13]," y.z^3 + ",self.actualPoly[14]," x.z^3",)

    def printAppxPoly(self):
        print(self.appxPoly[0]," x^2 + ",self.appxPoly[1]," y^2 + ",self.appxPoly[2]," z^2 + ",self.appxPoly[3]," x.y ",self.appxPoly[4]," yz ",self.appxPoly[5]," x.z ")

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

        P.printAppxPoly()
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
        #polyAppx=(P.appxPoly[0]*(x*x))+(P.appxPoly[1]*(y*y))+(P.appxPoly[2]*(z*z))+(P.appxPoly[3]*(x*y))+(P.appxPoly[4]*(y*z))+(P.appxPoly[5]*(x*z))
        polyAppx=(P.appxPoly[0]*(x*x))
        #model.optimize()
        #print(polyAppx)
        model.addConstr(polyAppx>=7)

        model.setObjective(polyAppx)
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






p1=Polynomial([1,1,1,1,1.1,1,1,0.8,1,0.2,1,1,0.1,0.01,1],[[1,2],[1,2],[1,2]],4)
p2=Polynomial([0,0,0.1,0,0.2,0,0,0.9,1,0,1,0,0.1,0.01,0],[[-1,1],[0,1],[-1,0]],4)
print("")
v2=Verification(p2,3000)
startTime=time.time()
v2.isSafe()
print("\n\n=======================================")
print("Result: ",v2.isSafe())
print("Time Taken: ",time.time()-startTime)
#p2.polyPower(45).printAppxPoly()
