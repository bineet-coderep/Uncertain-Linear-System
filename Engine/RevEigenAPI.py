'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System, Perturbation in the Eigenvectors and Eigenvalues,
compute the uncertainty that is introduced in the given Linear System.

Documentation: Not yet available. (TODO)
'''

import numpy as np
import numpy.linalg as LA
from IntervalNormAPI import *

PACE='fast'


class RevEigenDecomp:
    '''
    Given a Linear System A, Eigenvector perturbation dEigVal, Eigenvalue
    perturbation dEigVec, compute the uncertainty in A.
    '''

    def __init__(self,A,dEVal,dEVec):
        self.A=A # Linear System, numpy array.
        self.dEigVal=dEVal # Dictionary, Perturbation in Eigenvalues
        self.dEigVec=dEVec # Dictionary, Perturbation in Eigenvectors
        '''
        Following dictionary data-structure has been used to represent
        a error matrix:
        {
            (i,j): [a,b]
        }
        Indicating array[i][j] can pick any value within the range [a,b]
        '''
        self.n=A.shape[0] # Dimension of the System


    @staticmethod
    def matrixify(mat,n):
        Ac=np.zeros((n,n),dtype=object)
        for i in range(n):
            for j in range(n):
                Ac[i][j]=mpi(0,0)
        for key in mat:
            Ac[key[0]][key[1]]=mpi(mat[key][0],mat[key][1])
        return Ac

    @staticmethod
    def matrixify2(mat,n):
        Ac=np.zeros(n,dtype=object)
        for i in range(n):
            Ac[i]=mpi(0,0)
        for key in mat:
            Ac[key]=mpi(mat[key][0],mat[key][1])
        return Ac

    @staticmethod
    def dictionarify(mat):
        Er={}
        n=mat.shape[0]
        for i in range(n):
            for j in range(n):
                a=float(nstr(mat[i][j]).split(',')[0][1:])
                b=float(nstr(mat[i][j]).split(',')[1][:-1])
                if (a!=0 and b!=0):
                    Er[(i,j)]=[a,b]
        return Er


    @staticmethod
    def getPertEigValues(r,d):
        n=r.shape[0]
        q=r+RevEigenDecomp.matrixify2(d,n)
        Q=np.zeros((r.shape[0],r.shape[0]),dtype=object)
        for i in range(n):
            Q[i][i]=q[i]
        return Q


    def getUncertainty(self):
        '''
        Returns the uncertainty introduced in A due to dEigVal and dEigVec.
        Returns: (Norm of the Error Matrix, Interval Error Matrix)
        '''
        print("Started: Eigen Decomposition.....")
        (eigVal,eigVec)=LA.eig(self.A)
        print(".....Ended: Eigen Decomposition")
        pEVal=RevEigenDecomp.getPertEigValues(eigVal,self.dEigVal)
        pEVec=(eigVec.real+RevEigenDecomp.matrixify(self.dEigVec,self.n))
        print("Started: Inverting.....")
        pA=np.matmul(np.matmul(pEVec,pEVal),IntervalMatrix(pEVec).inverse())
        print(".....Started: Inverting")
        Er=pA-self.A
        print("Started: Norm.....")
        normEr=IntervalNorm(RevEigenDecomp.dictionarify(Er),self.n,PACE).getNorm()
        return (normEr,Er)
        print(".....Ended: Norm")

    def getReport(self):
        (normE,E)=self.getUncertainty()
        #norm_dEVal=IntervalNorm(self.dEigVal,self.n,PACE).getNorm()
        norm_dEVal=RevEigenDecomp.eigValNorm(self.dEigVal,self.n)
        print("Started: Norm.....")
        norm_dEVec=IntervalNorm(self.dEigVec,self.n,PACE).getNorm()
        print(".....Ended: Norm")
        return (norm_dEVal,norm_dEVec,normE,E)

    @staticmethod
    def eigValNorm(dVal,n):
        s=0
        for key in dVal:
            s=dVal[key][1]-dVal[key][0]
        s=s/n
        return s

    def printReport(self):
        (norm_dEVal,norm_dEVec,normE,E)=self.getReport()
        print("Started: Norm.....")
        if PACE=='fast':
            normA=LA.norm(self.A,ord='fro')
        else:
            normA=LA.norm(self.A,ord=2)
        print(".....Ended: Norm")
        print("-------Report-------")
        print("Norm of A: ",normA)
        print("Norm of E: ",normE)
        print("Percentage of Error (Norm): ",(normE*100)/normA)
        print("Average Perturbation in Eigenvalue: ",norm_dEVal)
        print("Norm of Perturbation in Eigenvector: ",norm_dEVec)
        print("Error Matrix E: ")
        print(E)
        print("------End of Report------")

# Tester --------------

if False:
    A=np.array([
    [1,0,0,0,0],
    [0,1,0,0,0],
    [0,0,1,0,0],
    [0,0,0,1,0],
    [0,0,0,0,1],
    ])
    dEVal={0: [-0.1,0.1]}
    dEVec={(3,0): [-0.1,0.1]}
    eig=RevEigenDecomp(A,dEVal,dEVec)
    eig.printReport()
