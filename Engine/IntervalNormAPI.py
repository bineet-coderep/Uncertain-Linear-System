'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

Most the functions are based on the paper:
- 'Norms of Interval Matrices'
by Raena Farhadsefat, Jirı Rohn and Taher Lotf


Documentation: Not yet available. (TODO)
'''
from mpmath import *
import numpy.linalg as LA
import numpy as np

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
        #print("n: ",n)
        #print(mat)
        ret=np.zeros((n,n),dtype=object)
        for i in range(n):
            for j in range(n):
                if (isinstance(mat[i][j],int)) or (isinstance(mat[i][j],float)):
                    ret[i][j]=(mat[i][j],mat[i][j])
                else:
                    #print("i: ",mat[i][j],type(mat[i][j]))
                    ret[i][j]=(float(nstr(mat[i][j]).split(',')[0][1:]),float(nstr(mat[i][j]).split(',')[1][:-1]))
        #print("Ret: ",ret)
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
            return np.asarray([[mpi(IntervalMatrix.intervalDiv(m[1][1],determinant)), mpi(IntervalMatrix.intervalDiv((-1*m[0][1][0],-1*m[0][1][1]),determinant))],
                    [mpi(IntervalMatrix.intervalDiv((-1*m[1][0][0],-1*m[1][0][1]),determinant)), mpi(IntervalMatrix.intervalDiv(m[0][0],determinant))]])

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

class IntervalNorm:
    '''
    Computes Norm-2 of a Interval Matrix
    based on the paper 'Norms of Interval Matrices'
    by Raena Farhadsefat, Jirı Rohn and Taher Lotf
    '''

    def __init__(self,matrix,s,p='slow'):
        self.A=matrix
        self.pace=p
        self.n=s

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
        Break the error matrix to A=[Ac-delta,Ac+delta)
        '''
        Ac=np.zeros((self.n,self.n))
        delta=np.zeros((self.n,self.n))
        for key in self.A:
            Ac[key[0]][key[1]]=(self.A[key][0]+self.A[key][1])/2
            delta[key[0]][key[1]]=self.A[key][1]-Ac[key[0]][key[1]]
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
            y=IntervalNorm.generateSignBits(i,self.n,1)
            for j in range(2**self.n):
                z=IntervalNorm.generateSignBits(j,self.n,0)
                tmp=IntervalNorm.spectralNorm(Ac+(np.matmul(y,z)*delta))
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
