#from Unsafe import *
#from Visualize import *

import pickle

import sys,os

import numpy as np
import math

from Parameters import *
from SplitMet import *


class BicycleModel:
    '''
    Obtain the robot and cloud model reachable set
    '''

    @staticmethod
    def createMatrix(A,B,mode,h):
        ''' Creates a single matrix based on
        . or +.
        In case of . a rough approximation is
        done'''

        n1=np.size(A,0)
        if (np.size(B)>0):
            n2=np.size(B,1)
        else:
            n2=0
        n=n1+n2

        C=np.zeros((n,n),dtype=np.float)
        if mode=='+':
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1
        elif mode=='.':
            I=np.zeros((n1,n1),dtype=np.float)
            for i in range(n1):
                I[i][i]=1
            A2=h*A
            A2=np.add(I,A2)
            B2=h*B
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A2[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B2[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C

    def getDynamics(v, phi, delta):
        NX = 4  # x = x, y, v, yaw
        NU = 2  # a = [accel, steer]
        DT = 0.2  # [s] time tick
        WB = 2.5  # [m]

        if v==0:
            v=0.000001

        A = np.zeros((NX, NX))
        A[0, 0] = 1.0
        A[1, 1] = 1.0
        A[2, 2] = 1.0
        A[3, 3] = 1.0
        A[0, 2] = DT * math.cos(phi)
        A[0, 3] = - DT * v * math.sin(phi)
        A[1, 2] = DT * math.sin(phi)
        A[1, 3] = DT * v * math.cos(phi)
        A[3, 2] = DT * math.tan(delta) / WB

        B = np.zeros((NX+1, NU+1))
        B[0, 2]= DT * v * math.sin(phi) * phi # For C[0]
        B[1, 2]= - DT * v * math.cos(phi) * phi # For C[1]
        B[2, 0] = DT
        B[3, 1] = DT * v / (WB * math.cos(delta) ** 2)
        B[3, 2] = - DT * v * delta / (WB * math.cos(delta) ** 2) # For C[3]
        B[4, 2] = 1

        p=220
        Er={
        #(0,0): [1-(p/100),1+(p/100)],
        (0,2): [1-(p/100),1+(p/100)],
        (0,3): [1-(p/100),1+(p/100)],
        (0,6): [1-(p/100),1+(p/100)],
        #(1,1): [1-(p/100),1+(p/100)],
        (1,2): [1-(p/100),1+(p/100)],
        (1,3): [1-(p/100),1+(p/100)],
        (1,6): [1-(p/100),1+(p/100)],
        }
        #print(A[0, 2], A[0, 3], B[0, 2])

        return A, B, Er

    def getReachSets(X, Y, yawList, Vel, D, A):
        '''
        Pickle a set of Robot-Cloud pair, given a initial set,
        and U
        '''
        rcPair={}
        rbt_all=[]
        cld_all=[]
        dim=7
        C=[0]*dim
        V=np.identity(dim)

        rg=int(math.floor(len(D)/((len(D)*PER_COVERAGE)/100)))
        listT=list(range(1,len(D),rg))
        listT.append(len(D)-1)

        rsList=[]
        for i in range(len(D)):
            if (i in listT):
                print(i)
                P=[(1,1)]*dim
                P[0]=(X[i],X[i]) # Update x
                P[1]=(Y[i],Y[i]) # Update y
                P[2]=(Vel[i],Vel[i]) # Update v
                P[3]=(yawList[i],yawList[i]) # Update \phi
                P[4]=(A[i],A[i]) # Update a
                P[5]=(D[i],D[i]) # Update \delta
                initialSet=(C,V,P)
                reachS=BicycleModel.getRSPoint(initialSet)
                rsList.append(reachS)

        return rsList

    def getRSPoint(initialSet):
        C=initialSet[0]
        V=initialSet[1]
        P=initialSet[2]
        v=P[2][0]
        phi=P[3][0]
        delta=P[5][0]
        (dynA,dynB,Er)=BicycleModel.getDynamics(v,phi,delta)
        A=BicycleModel.createMatrix(dynA,dynB,'+',1)
        rs=Split(A,Er,initialSet,1)
        reachS=rs.getReachableSet()
        return reachS
