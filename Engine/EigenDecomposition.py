'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

IN PROGRESS

- Guivers Method.
- Explored scalability of interval arithmetic based method

Future Notes:
Most the functions are based on the paper:
- 'A note on the eigenvectors of perturbed matrices with applications to linear positive system'
by Guiver et. al

- Linear Dynamical System: dot{x} = (A+BPC)x; where P (BPC) is the perturbation.
- Given a perturbation, find out the eigenvalues and eigenvectors of the perturbed matrix.

Documentation: Not yet available. (TODO)
'''


from mpmath import *
import numpy as np
import numpy.linalg as LA
import time




class IntervalMatrix:
    '''
    Operations on interval matrices
    '''

    def __init__(self,mtrx,p='slow'):
        self.mat=mtrx
        '''
        A numpy array with elements as mpi(a,b)
        '''
        self.pace=p
        self.n=mtrx.shape[0]

    @staticmethod
    def formatize(mat):
        n=mat.shape[0]
        ret=np.zeros((n,n),dtype=object)
        for i in range(n):
            for j in range(n):
                if type(mat[i][j]) is int:
                    ret[i][j]=(mat[i][j],mat[i][j])
                else:
                    ret[i][j]=(float(nstr(mat[i][j]).split(',')[0][1:]),float(nstr(mat[i][j]).split(',')[1][:-1]))
        return ret

    @staticmethod
    def intervalAdd(i1,i2):
        return (i1[0]+i2[0],i1[1]+i2[1])

    @staticmethod
    def intervalSub(i1,i2):
        return (i1[0]-i2[1],i1[1]-i2[0])

    @staticmethod
    def intervalDiv(i1,i2):
        return IntervalMatrix.intervalMul(i1,(1/i2[0],1/i2[1]))

    @staticmethod
    def intervalMul(i1,i2):
        return (min(i1[0]*i2[0],i1[0]*i2[1],i1[1]*i2[0],i1[1]*i2[1]),max(i1[0]*i2[0],i1[0]*i2[1],i1[1]*i2[0],i1[1]*i2[1]))

    @staticmethod
    def display(m):
        r=m.shape[0]
        c=m.shape[1]
        for i in range(r):
            for j in range(c):
                print(m[i][j],"          ",end="")
            print()

    def inverse(self):
        '''
        Compute inverse of the matrix.
        Return: A numpy array with elements mpi(a,b)
        '''
        m=IntervalMatrix.formatize(self.mat)
        m=m.tolist()

        determinant = IntervalMatrix.getMatrixDeternminant(m)
        '''print(determinant)
        exit()'''
        #special case for 2x2 matrix:
        if len(m) == 2:
            return np.asarray([[IntervalMatrix.intervalDiv(m[1][1],determinant), IntervalMatrix.intervalDiv((-1*m[0][1][0],-1*m[0][1][1]),determinant)],
                    [IntervalMatrix.intervalDiv((-1*m[1][0][0],-1*m[1][0][1]),determinant), IntervalMatrix.intervalDiv(m[0][0],determinant)]])

        #find matrix of cofactors
        cofactors = []
        for r in range(len(m)):
            cofactorRow = []
            for c in range(len(m)):
                minor = IntervalMatrix.getMatrixMinor(m,r,c)
                cofactorRow.append((((-1)**(r+c)) * IntervalMatrix.getMatrixDeternminant(minor)[0],((-1)**(r+c)) * IntervalMatrix.getMatrixDeternminant(minor)[1]))
            cofactors.append(cofactorRow)

        n=len(m)
        cofactorsTranspose=[]
        for i in range(n):
            l=[]
            for j in range(n):
                l.append((cofactors[j][i][0],cofactors[j][i][1]))
            cofactorsTranspose.append(l)

        inverseList=[]

        for r in range(len(cofactors)):
            l=[]
            for c in range(len(cofactors)):
                l.append(IntervalMatrix.intervalDiv(cofactorsTranspose[r][c],determinant))
            inverseList.append(l)

        inv=np.zeros((n,n),dtype=object)
        for i in range(n):
            for j in range(n):
                inv[i][j]=mpi(inverseList[i][j][0],inverseList[i][j][1])

        return inv

    @staticmethod
    def getMatrixMinor(m,i,j):
        return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

    @staticmethod
    def getMatrixDeternminant(m):
        #base case for 2x2 matrix
        if len(m) == 2:
            return IntervalMatrix.intervalSub(IntervalMatrix.intervalMul(m[0][0],m[1][1]),IntervalMatrix.intervalMul(m[0][1],m[1][0]))

        determinant = (0,0)
        for c in range(len(m)):
            determinant = IntervalMatrix.intervalAdd((IntervalMatrix.intervalMul((((-1)**c)*m[0][c][0],((-1)**c)*m[0][c][1]),IntervalMatrix.getMatrixDeternminant(IntervalMatrix.getMatrixMinor(m,0,c)))),determinant)
        return determinant

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
        Break the matrix to mat=[Ac-delta,Ac+delta)
        '''
        m=IntervalMatrix.formatize(self.mat)
        Ac=np.zeros((self.n,self.n))
        delta=np.zeros((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                Ac[i][j]=(m[i][j][0]+m[i][j][1])/2
                delta[i][j]=m[i][j][1]-Ac[i][j]
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
            y=IntervalMatrix.generateSignBits(i,self.n,1)
            for j in range(2**self.n):
                z=IntervalMatrix.generateSignBits(j,self.n,0)
                tmp=IntervalMatrix.spectralNorm(Ac+(np.matmul(y,z)*delta))
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


class EigenDecompose:
    '''
    This class computes the eigenvalues of the perturbed matrix using
    Bauer-Fike Theorem. And the eigenvectors are computed based on the
    techniques in the paper "A note on the eigenvectors of perturbed
    matrices with applications to linear positive system" by Guiver et. al.
    '''

    def __init__(self,A,B,P,C):
        '''
        A_pert=A+BPC
        Please refer to corrollary 2.4 for details
        '''
        self.A=A
        self.b=B
        self.C=C
        self.q=P

    def getEigenValues(self,pace='slow'):
        '''
        The Eigenvalues are obtained using Bauer-Fike Theorem
        '''
        E=np.matmul(np.matmul(self.b,self.q),self.C) # Error Matrix
        (A_EV,A_EVec)=LA.eig(self.A) # Compute the Eigenvalues and Eigenvectors of A
        k2=EigenDecompose.conditionNum(A_EVec) # Computes the Condition Number
        normE=IntervalMatrix(E,pace).getNorm() # Computes Interval Norm of interval matrix E
        distance=k2*normE # Computes the distance according to Bauer-Fike Thorem
        EigenVals=np.zeros(A_EV.shape[0],dtype=object)
        for i in range(A_EV.shape[0]):
            EigenVals[i]=mpi(A_EV[i]-distance,A_EV[i]+distance)
        return EigenVals

    def decompose(self):
        '''
        The Eigenvectors are obtained using corrollary 2.4(b) of the
        above mentioned paper
        '''

        evals=self.getEigenValues()
        n=self.A.shape[0]
        evects=np.zeros((n,n),dtype=object)
        i=0
        for e in evals:
            lI=EigenDecompose.getLamdaI(e,n) # e \time I Interval matrix
            t=IntervalMatrix(lI-self.A).inverse() # (\lambda*I - A)^-1
            qC=np.matmul(self.q,self.C)
            evec=np.matmul(qC,t)
            evects[i]=evec
            i=i+1
        evects=evects.transpose()
        return (evals,evects)

    @staticmethod
    def getLamdaI(e,n):
        lI=np.zeros((n,n),dtype=object)
        for i in range(n):
            lI[i]=e
        return lI

    @staticmethod
    def conditionNum(V):
        return LA.norm(V,ord=2)*LA.norm(LA.inv(V),ord=2)



# Tester ------------------------
if False:
    A=np.array([
    [3,0,0,0,0,2,4],
    [1,2,3,0,0,2.9,2.9],
    [8,1,2,0,0,2.9,2.9],
    [7,0,0,8,2,3.9,3.9],
    [8,0,0,3,7,3.9,3.9],
    [0,0,0,0,0,6,3],
    [0,0,0,0,0,2,1],
    ])
    b=np.array([
    [0],
    [0],
    [1],
    [0],
    [0],
    [1],
    [0]
    ])
    q=np.array([[0,0,mpi(2,3),0,0,mpi(2,3),0]])
    C=np.array([
    [1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0],
    [0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0],
    [0,0,0,1,0,0,0],
    [0,1,0,0,0,0,0],
    [0,0,0,0,1,0,0],
    ])

    d=EigenDecompose(A,b,q,C)
    (v,vect)=d.decompose()
    print(v)
    print("")
    print(vect)
