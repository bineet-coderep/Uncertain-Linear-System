'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- POC on  Guivers Method.
- Explored scalability of interval arithmetic based method with mpmath
'''

from mpmath import *
import numpy as np
import time

mp.dps = 15

class IntervalMatrix:
    def __init__(self,mtrx):
        self.mat=mtrx
        '''
        A numpy array with elements as mpi(a,b)
        '''

    def inverse(self):
            '''
            Compute inverse of the matrix.
            Return: A numpy array with elements mpi(a,b)
            '''
            m=self.mat.tolist()

            determinant = IntervalMatrix.getMatrixDeternminant(m)
            #special case for 2x2 matrix:
            if len(m) == 2:
                return np.asarray([[m[1][1]/determinant, -1*m[0][1]/determinant],
                        [-1*m[1][0]/determinant, m[0][0]/determinant]])

            #find matrix of cofactors
            cofactors = []
            for r in range(len(m)):
                cofactorRow = []
                for c in range(len(m)):
                    minor = IntervalMatrix.getMatrixMinor(m,r,c)
                    cofactorRow.append(((-1)**(r+c)) * IntervalMatrix.getMatrixDeternminant(minor))
                cofactors.append(cofactorRow)
            cofactors = np.transpose(np.asarray(cofactors)).tolist()
            for r in range(len(cofactors)):
                for c in range(len(cofactors)):
                    cofactors[r][c] = cofactors[r][c]/determinant


            inv=np.asarray(cofactors)
            return inv


    @staticmethod
    def getMatrixMinor(m,i,j):
        return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

    @staticmethod
    def getMatrixDeternminant(m):
        #base case for 2x2 matrix
        if len(m) == 2:
            return m[0][0]*m[1][1]-m[0][1]*m[1][0]

        determinant = 0
        for c in range(len(m)):
            determinant += ((-1)**c)*m[0][c]*IntervalMatrix.getMatrixDeternminant(IntervalMatrix.getMatrixMinor(m,0,c))
        return determinant


ABI=np.array([
[mpi(2,2.01),2,1,0,2,3,0,1,0,3],
[1,mpi(2,2.01),1,1,0,3,0,-1,0,3],
[0,2,mpi(2,2.01),0,2,-2,0,1,0,3],
[-3,2,1,mpi(2,2.01),0,3,0,1,0,0],
[2,0,1,2,mpi(2,2.01),3,0,4,0,3],
[-2,1,1,0,-2,mpi(2,2.01),0,-1,3,0],
[1,2,1,0,2,-1,mpi(2,2.01),1,0,3],
[0,2,1,0,2,3,0,mpi(2,2.01),0,3],
[-1,1,1,0,2,3,0,1,mpi(2,2.01),3],
[-2,2,1,0,2,3,0,1,0,mpi(2,2.01)]
])

ABU=np.array([
[2,2,1,0,2,3,0,1,0,3],
[1,2,1,1,0,3,0,-1,0,3],
[0,2,2,0,2,-2,0,1,0,3],
[-3,2,1,2,0,3,0,1,0,0],
[2,0,1,2,2,3,0,4,0,3],
[-2,1,1,0,-2,2,0,-1,3,0],
[1,2,1,0,2,-1,2,1,0,3],
[0,2,1,0,2,3,0,2,0,3],
[-1,1,1,0,2,3,0,1,2,3],
[-2,2,1,0,2,3,0,1,0,2]
])

M1I=np.array([
[mpi(1,1),6,3,4],
[mpi(2,2),mpi(0,0),1,2],
[mpi(1,1),0,mpi(0,0),1],
[mpi(3,3),0,3,mpi(0,0)],
])

M1=np.array([
[1,6,3,4],
[2,0,1,2],
[1,0,0,1],
[3,0,3,0],
])

print(np.linalg.inv(M1))
print("")
matA=IntervalMatrix(M1)
start_time=time.time()
invA=matA.inverse()
end_time=time.time()-start_time
print(invA)
print(end_time)
