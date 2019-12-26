'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

IN PROGRESS

- POC on  Guivers Method.
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
import time

mp.dps = 50


class IntervalMatrix:
    '''
    Operations on interval matrices
    '''

    def __init__(self,mtrx):
        self.mat=mtrx
        '''
        A numpy array with elements as mpi(a,b)
        '''


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
        m=self.mat.tolist()

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
                inv[i][j]=(inverseList[i][j][0],inverseList[i][j][1])

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


#A=np.array([[2,2],[2,mpi(3,4)]])
AB=np.array([
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

ABI=np.array([
[(2,2),(2,2),(1,1),(0,0),(2,2),(3,3),(0,0),(1,1),(0,0),(3,3)],
[(1,1),(2,2),(1,1),(1,1),(0,0),(3,0),(0,0),(-1,-1),(0,0),(3,3)],
[(0,0),(2,2),(2,2),(0,0),(2,2),(-2,-2),(0,0),(1,1),(0,0),(3,3)],
[(-3,-3),(2,2),(1,1),(2,2),(0,0),(3,3),(0,0),(1,1),(0,0),(0,0)],
[(2,2),(0,0),(1,1),(2,2),(2,2),(3,3),(0,0),(4,4),(0,0),(3,3)],
[(-2,-2),(1,1),(1,1),(0,0),(-2,-2),(2,2),(0,0),(-1,-1),(3,3),(0,0)],
[(1,1),(2,2),(1,1),(0,0),(2,2),(-1,-1),(2,2),(1,1),(0,0),(3,3)],
[(0,0),(2,2),(1,1),(0,0),(2,2),(3,3),(0,0),(2,2),(0,0),(3,3)],
[(-1,-1),(1,1),(1,1),(0,0),(2,2),(3,3),(0,0),(1,1),(2,2),(3,3)],
[(-2,-2),(2,2),(1,1),(0,0),(2,2),(3,3),(0,0),(1,0),(0,0),(2,2)]
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
AC=np.array([
[1,2,1,0,2,3,0,1,0,3],
[1,1,1,1,0,3,0,-1,0,3],
[0,2,1,0,2,-2,0,1,0,3],
[-3,2,1,1,0,3,0,1,0,0],
[2,0,1,2,1,3,0,4,0,3],
[-2,1,1,0,-2,1,0,-1,3,0],
[1,2,1,0,2,-1,1,1,0,3],
[0,2,1,0,2,3,0,1,0,3],
[-1,1,1,0,2,3,0,1,1,3],
[-2,2,1,0,2,3,0,1,0,1]
])
A=np.array([
[1,2,1],
[1,1,1],
[0,2,1]
])
B=np.array([
[(2,2),(2,2),(1,1)],
[(1,1),(2,2),(1,1)],
[(0,0),(2,2),(2,2)]
])
BC=np.array([
[2,2,1],
[1,2,1],
[0,2,2]
])
'''B=np.array([
[(2,3),(2,2)],
[(1,1),(2,3)],
])'''

M1=np.array([
[1,6,3,4],
[2,0,1,2],
[1,0,0,1],
[3,0,3,0],
])

M1I=np.array([
[(1,1),(6,6),(3,3),(4,4)],
[(2,2),(0,0),(1,1),(2,2)],
[(1,1),(0,0),(0,0),(1,1)],
[(3,3),(0,0),(3,3),(0,0)],
])










#--------------------------------------------
#print(np.linalg.inv(M1))

print("")
matA=IntervalMatrix(ABI)
start_time=time.time()
invA=matA.inverse()
end_time=time.time()-start_time
print(invA)
print(end_time)

#print("-")
#IntervalMatrix.display(invA)
#matA.inverse()
