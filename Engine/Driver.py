'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

from VerificationEngine import *
from Benchmarks import *
from VisualizeAPI import *
from RobustMetricAPI import *

class DriverBloat:


    @staticmethod
    def createMatrix(A,B,mode,h=0):
        ''' Creates a single matrix based on
        . or +.
        In case of . a roungh approximation is
        done'''

        n1=np.size(A,0)
        if (np.size(B)>0):
            n2=np.size(B,1)
        else:
            n2=0
        n=n1+n2

        C=np.zeros((n,n))
        for i in range(n1):
            for j in range(n1):
                C[i][j]=A[i][j]
        for i in range(n1):
            j2=0
            for j in range(n1,n1+n2):
                C[i][j]=B[i][j2]
                j2=j2+1

        if mode=='+':
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C

    @staticmethod
    def bpcP2E(A,B,mode,h,b,p,c):
        '''
        Given error matrices B, P, C in percent error,
        converts it into E dictionary with absolute value.

        [0.9,1.1] means +-10% ([90%,110%]) perturbation.
        '''
        E={}
        conA=DriverBloat.createMatrix(A,B,mode,h)
        n=conA.shape[0]
        ErP=np.matmul(np.matmul(b,p),c)

        for i in range(n):
            for j in range(n):
                a=float(nstr(ErP[i][j]).split(',')[0][1:])
                b=float(nstr(ErP[i][j]).split(',')[1][:-1])
                if b>a:
                    i1=conA[i][j]*a
                    i2=conA[i][j]*b
                    E[(i,j)]=[i1,i2]

        return E

    @staticmethod
    def stableSystem1():
        A=Benchmarks.StableSystem1.A
        B=Benchmarks.StableSystem1.B
        mode='.'
        E={
        (0,2): [-0.1,0.1],
        (2,1): [-1,1],
        }
        IS_C=np.array([0,0,0])
        IS_V=np.array([
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        #vrfy.plotTime(0,3,0.01)
        vrfy.plotTimeCompare(0,30,1,['Kagstrom1','Kagstrom2','Loan'],'slow')

    def stableSystem2():
        A=Benchmarks.StableSystem2.A
        B=Benchmarks.StableSystem2.B
        mode='.'
        '''E={
        (0,1): [-0.1,0.1]
        }'''
        b=np.array([
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0],
        [1,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0])
        IS_V=np.array([
        [0,1,0],
        [0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5]
        method='Kagstrom2'


        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        #vrfy.plotTime(0,3,0.01)
        vrfy.plotTimeCompare(0,t,1,['Kagstrom2','Kagstrom1','Loan'],'slow')

    def stableSystem3():
        A=Benchmarks.StableSystem3.A
        B=Benchmarks.StableSystem3.B
        mode='.'
        '''E={
        (0,1): [-0.1,0.1]
        }'''
        b=np.array([
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0],
        [1,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0])
        IS_V=np.array([
        [0,1,0],
        [0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        #vrfy.plotTime(0,3,0.01)
        vrfy.plotTimeCompare(0,30,1,['Kagstrom2','Kagstrom1','Loan'],'slow')

    def stableSystem4():
        A=Benchmarks.StableSystem4.A
        B=Benchmarks.StableSystem4.B
        mode='.'
        E={
        (0,2): [-0.1,0.1],
        (2,1): [-1,1],
        }
        IS_C=np.array([0,0,0])
        IS_V=np.array([
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        #vrfy.plotTime(0,3,0.01)
        vrfy.plotTimeCompare(0,30,0.1,['Kagstrom2','Kagstrom1','Loan'],'slow')

    @staticmethod
    def illustExample():
        A=Benchmarks.IllustExample.A
        B=Benchmarks.IllustExample.B
        mode='.'
        '''E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0],
        [1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,1,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0,0,0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,1]
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5,5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        #vrfy.plotTime(0,3,0.01)
        vrfy.plotTimeCompare(0,20,1,['Kagstrom1','Kagstrom2','Loan'],'fast')

    @staticmethod
    def flightEnvelope():
        A=Benchmarks.FlightEnvelope.A
        B=Benchmarks.FlightEnvelope.B
        mode='.'
        '''E={
        (3,7): [-1,1],
        (4,6): [-1,1],
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [1],
        [0],
        [0],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,0,0,mpi(-0.1,0.1),0,0,0,0,0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        #vrfy.plotTime(0,10,0.01,'fast')
        vrfy.plotTimeCompare(0,20,1,['Kagstrom1','Kagstrom2','Loan'],'fast')

    @staticmethod
    def coOPVehiclesI():
        A=Benchmarks.CoOPVehiclesI.A
        B=Benchmarks.CoOPVehiclesI.B
        mode='.'
        '''E={
        (2,0): [-0.01,0.01],
        (2,1): [-0.1,0.1],
        (2,8): [-0.001,0.001],
        (5,0): [-0.01,0.01],
        (5,8): [-0.01,0.01]
        }
        linE={
        (0,9): [-0.01,0.01],
        (1,9): [-0.1,0.1],
        (2,9): [-0.001,0.001],
        (3,9): [-0.01,0.01],
        (4,9): [-0.01,0.01]
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [1,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0,0,0,0,0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,1]
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5,5,5,5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        vrfy.plotTime(0,20,1)
        vrfy.plotTimeCompare(0,20,1,['Kagstrom1','Kagstrom2','Loan'])

    def coOPVehiclesII():
        A=Benchmarks.CoOPVehiclesII.A
        B=Benchmarks.CoOPVehiclesII.B
        mode='.'
        '''E={
        (2,0): [-0.01,0.01],
        (2,1): [-0.1,0.1],
        (2,8): [-0.001,0.001],
        (5,0): [-0.01,0.01],
        (5,8): [-0.01,0.01]
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [1,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0,0,0,0,0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,1]
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5,5,5,5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        vrfy.plotTimeCompare(0,20,1,['Kagstrom1','Kagstrom2','Loan'])

    @staticmethod
    def pkpd2():
        A=Benchmarks.PKPD2.A
        B=Benchmarks.PKPD2.B
        mode='.'
        '''linE={
        (0,4): [-0.01,0.01],
        }
        E={
        (3,3): [-0.01,0.01],
        }'''
        b=np.array([
        [0],
        [1],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0],
        [0,0,0,0,0],
        [0,0,1,0,0],
        [0,0,0,0,0],
        [0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0],
        [0,0,1,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,1,0],
        [0,0,0,0,0,1]
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5]
        method='Loan'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)

        #vrfy.plotTime(0,50,1)
        vrfy.plotTimeCompare(0,20,1,['Kagstrom1','Kagstrom2','Loan'])

    @staticmethod
    def dcConv():
        A=Benchmarks.DCConv.A
        B=Benchmarks.DCConv.B
        mode='.'
        E={
        (1,1): [-2,2],
        }
        IS_C=np.array([0,0,0])
        IS_V=np.array([
        [0,1,0,0,],
        [0,0,1,0,],
        [0,0,0,1,],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        vrfy.plotTimeCompare(0,50,0.001,['Kagstrom1','Loan'])
        #vrfy.plotTime(0,3,0.01)

    @staticmethod
    def spaceCraftRndzvs():
        A=Benchmarks.SpaceCraftRndzvs.A
        B=Benchmarks.SpaceCraftRndzvs.B
        mode='.'
        '''E={
        (2,0): [-1,1],
        (2,3): [-0.1,0.1],
        (3,2): [-0.5,0.5],
        }'''
        b=np.array([
        [1],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1)]])
        C=np.array([
        [1,0,0,0,0,0],
        [0,1,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,1,0],
        [0,1,0,0,0,0],
        [0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0,0],
        [0,0,1,0,0,0,0],
        [0,0,0,1,0,0,0],
        [0,0,0,0,1,0,0],
        [0,0,0,0,0,1,0],
        [0,0,0,0,0,0,1]
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        #vrfy.plotTime(0,20,1)
        vrfy.plotTimeCompare(0,20,1,['Kagstrom1','Kagstrom2','Loan'])

    @staticmethod
    def holesCXc():
        A=Benchmarks.HolesCXc.A
        B=Benchmarks.HolesCXc.B
        mode='.'
        '''E={
        (0,3): [-0.5,0.5],
        (1,2): [-0.1,0.1],
        (3,2): [-0.005,0.005],
        (4,3): [-0.0045,0.0045]
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [1,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS_C=np.array([0,0,0,0,0,0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5,5,5,5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        vrfy.plotTimeCompare(0,20,1,['Kagstrom1','Kagstrom2','Loan'])

    @staticmethod
    def holesPDp():
        A=Benchmarks.HolesPDp.A
        B=Benchmarks.HolesPDp.B
        mode='.'
        E={
        (0,1): [-0.001,0.001],
        (1,1): [-0.001,0.001],
        (1,0): [-0.1,0.1]
        }
        IS_C=np.array([0,0])
        IS_V=np.array([
        [0,1,0],
        [0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5]
        method='Loan'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        vrfy.plotTime(0,15,0.01)

    @staticmethod
    def holesPXp():
        A=Benchmarks.HolesPXp.A
        B=Benchmarks.HolesPXp.B
        mode='.'
        E={
        (0,0): [-0.1,0.1],
        (0,1): [-0.1,0.1],
        (0,2): [-0.1,0.1],
        (2,0): [-0.001,0.001],
        (2,2): [-0.001,0.001],
        (4,4): [-0.0001,0.0001]
        }
        IS_C=np.array([0,0,0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5,5,5]
        method='Loan'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        vrfy.plotTime(0,3,0.01)

    @staticmethod
    def motorTransmission1():
        A=Benchmarks.MotorTransmission1.A
        B=Benchmarks.MotorTransmission1.B
        mode='.'
        '''E={
        (0,6): [-0.1,0.1],
        (1,6): [-1,1]
        }'''
        b=np.array([
        [1],
        [0],
        [1],
        [0],
        [1],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0,0,0,0,0,0],
        [0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0],
        [0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)

        IS_C=np.array([0,0,0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5,5,5]
        method='Kagstrom2'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        #vrfy.plotTime(0,5,0.1)
        #vrfy.comparePace(['Kagstrom','Loan'])
        vrfy.plotTimeCompare(0,20,1,['Kagstrom1','Kagstrom2','Loan'])

    @staticmethod
    def motorTransmission2():
        A=Benchmarks.MotorTransmission1.A
        B=Benchmarks.MotorTransmission1.B
        mode='.'
        E={
        (5,0): [-0.1,0.1],
        (5,1): [-0.1,0.1]
        }
        IS_C=np.array([0,0,0,0,0])
        IS_V=np.array([
        [0,1,0,0,0,0],
        [0,0,1,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,1,0],
        [0,0,0,0,0,1],
        ])
        IS_r=2
        IS=[IS_C,IS_V,IS_r]
        t=20
        U=[5,5,5,5,5]
        method='loan'

        vrfy=VerifyBloat(A,B,E,IS,t,U,method)
        vrfy.plotTime(0,20,0.01)

class DriverDecomp:

    @staticmethod
    def bpcP2bpc(A,B,mode,h,b,p,c):
        '''
        Given error matrices B, P, C in percent error,
        converts it into E dictionary with absolute value.

        [0.9,1.1] means +-10% ([90%,110%]) perturbation.
        '''
        E={}
        conA=DriverBloat.createMatrix(A,B,mode,h)
        n=conA.shape[0]
        ErP=np.matmul(np.matmul(b,p),c)

        for i in range(n):
            for j in range(n):
                a=float(nstr(ErP[i][j]).split(',')[0][1:])
                b=float(nstr(ErP[i][j]).split(',')[1][:-1])
                if b>a:
                    i1=conA[i][j]*a
                    i2=conA[i][j]*b
                    E[(i,j)]=[i1,i2]

        return E

    @staticmethod
    def formatize(mat):
        n=mat.shape[0]
        ret=np.zeros((n,1),dtype=object)
        for i in range(n):
            #exit(0)
            if (isinstance(mat[i][0],int)) or (isinstance(mat[i][0],float)):
                ret[i][0]=(mat[i][0],mat[i][0])
            else:
                ret[i][0]=(float(nstr(mat[i][0]).split(',')[0][1:]),float(nstr(mat[i][0]).split(',')[1][:-1]))
        return ret

    @staticmethod
    def illustExample():
        A=Benchmarks.IllustExample.A
        B=Benchmarks.IllustExample.B
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0],
        [1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,1,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0]
        ])
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        ])

        print(np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        U=np.array([5,5,5,5,5,5,5,5])
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print(v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print(vP)

    @staticmethod
    def flightEnvelope():
        A=Benchmarks.FlightEnvelope.A
        B=Benchmarks.FlightEnvelope.B
        mode='.'
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [1],
        [0],
        [0],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,0,0,mpi(-0.1,0.1),0,0,0,0,0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]
        ])
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        print(np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print(v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print(vP)

    @staticmethod
    def stableSystem2():
        A=Benchmarks.StableSystem2.A
        B=Benchmarks.StableSystem2.B
        mode='.'
        '''E={
        (0,1): [-0.1,0.1]
        }'''
        b=np.array([
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0],
        [1,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        ])
        U=[]
        print(np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print(v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print(vP)

    @staticmethod
    def stableSystem3():
        A=Benchmarks.StableSystem3.A
        B=Benchmarks.StableSystem3.B
        mode='.'
        '''E={
        (0,1): [-0.1,0.1]
        }'''
        b=np.array([
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0],
        [1,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        ])
        U=[]
        print(np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print(v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print(vP)

    @staticmethod
    def coOPVehiclesII():
        A=Benchmarks.CoOPVehiclesII.A
        B=Benchmarks.CoOPVehiclesII.B
        mode='.'
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [1,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        ])
        U=[5,5,5,5,5,5,5,5,5,5]

        print(np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        print(vrfy.computeReachSetPertFree())
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print(v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print(vP)
        #vrfy.plotTimeCompare(0,10,0.01,['Kagstrom','Loan'])

    @staticmethod
    def pkpd2():
        A=Benchmarks.PKPD2.A
        B=Benchmarks.PKPD2.B
        U=[5,5,5,5,5]
        b=np.array([
        [0],
        [1],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0],
        [0,0,0,0,0],
        [0,0,1,0,0],
        [0,0,0,0,0],
        [0,0,0,0,0]
        ])
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1]
        ])

        print("Error: \n",np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print("Reach Set without Pert: \n",v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print("Reach Set with Pert: \n",vP)

    @staticmethod
    def spaceCraftRndzvs():
        A=Benchmarks.SpaceCraftRndzvs.A
        B=Benchmarks.SpaceCraftRndzvs.B
        mode='.'
        b=np.array([
        [1],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1)]])
        C=np.array([
        [1,0,0,0,0,0],
        [0,1,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,1,0],
        [0,1,0,0,0,0],
        [0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        print(np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print(v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print(vP)

    @staticmethod
    def holesCXc():
        A=Benchmarks.HolesCXc.A
        B=Benchmarks.HolesCXc.B
        mode='.'
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [1,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        print(np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print(v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print(vP)

    @staticmethod
    def motorTransmission1():
        A=Benchmarks.MotorTransmission1.A
        B=Benchmarks.MotorTransmission1.B
        mode='.'
        b=np.array([
        [1],
        [0],
        [1],
        [0],
        [1],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0,0,0,0,0,0],
        [0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0],
        [0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        print(np.matmul(np.matmul(b,q),C))
        print("\n\n\n\n")
        t=20
        vrfy=VerifyDecomp(A,B,b,C,q,IS,t,U)
        v=DriverDecomp.formatize(vrfy.computeReachSetPertFree())
        print(v)
        print("\n\n\n\n")
        vP=DriverDecomp.formatize(vrfy.computeReachSet())
        print(vP)

class DriverInterval:

    @staticmethod
    def bpcToE(b,q,C):
        Er=np.matmul(np.matmul(b,q),C)
        n=Er.shape[0]
        E={}
        for i in range(n):
            for j in range(n):
                a=float(nstr(Er[i][j]).split(',')[0][1:])
                b=float(nstr(Er[i][j]).split(',')[1][:-1])
                if (a!=0 or b!=0):
                    E[(i,j)]=[a,b]
        return E

    @staticmethod
    def illustExample():
        A=Benchmarks.IllustExample.A
        B=Benchmarks.IllustExample.B
        mode='.'
        '''E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0],
        [1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,1,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        t=20
        U=[5,5,5,5,5,5,5,5]
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

    @staticmethod
    def stableSystem2():
        A=Benchmarks.StableSystem2.A
        B=Benchmarks.StableSystem2.B
        mode='.'
        '''E={
        (0,1): [-0.1,0.1]
        }'''
        b=np.array([
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0],
        [1,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1]
        ])
        t=20
        U=[5,5]
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

    @staticmethod
    def stableSystem3():
        A=Benchmarks.StableSystem3.A
        B=Benchmarks.StableSystem3.B
        mode='.'
        '''E={
        (0,1): [-0.1,0.1]
        }'''
        b=np.array([
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0],
        [1,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1]
        ])
        t=20
        U=[5,5]
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

    @staticmethod
    def flightEnvelope():
        A=Benchmarks.FlightEnvelope.A
        B=Benchmarks.FlightEnvelope.B
        mode='.'
        '''E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [1],
        [0],
        [0],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,0,0,mpi(-0.1,0.1),0,0,0,0,0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]
        ])
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        E=DriverInterval.bpcToE(b,q,C)
        t=20
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

    @staticmethod
    def coOPVehiclesII():
        A=Benchmarks.CoOPVehiclesII.A
        B=Benchmarks.CoOPVehiclesII.B
        mode='.'
        '''E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [1,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        E=DriverInterval.bpcToE(b,q,C)
        t=20
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

    @staticmethod
    def pkpd2():
        A=Benchmarks.PKPD2.A
        B=Benchmarks.PKPD2.B
        mode='.'
        '''E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }'''
        b=np.array([
        [0],
        [1],
        [1],
        [0],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0],
        [0,0,0,0,0],
        [0,0,1,0,0],
        [0,0,0,0,0],
        [0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        E=DriverInterval.bpcToE(b,q,C)
        t=20
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

    @staticmethod
    def spaceCraftRndzvs():
        A=Benchmarks.SpaceCraftRndzvs.A
        B=Benchmarks.SpaceCraftRndzvs.B
        mode='.'
        '''E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }'''
        b=np.array([
        [1],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1)]])
        C=np.array([
        [1,0,0,0,0,0],
        [0,1,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,1,0],
        [0,1,0,0,0,0],
        [0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        E=DriverInterval.bpcToE(b,q,C)
        t=20
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

    @staticmethod
    def holesCXc():
        A=Benchmarks.HolesCXc.A
        B=Benchmarks.HolesCXc.B
        mode='.'
        '''E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }'''
        b=np.array([
        [0],
        [0],
        [1],
        [0],
        [0],
        [0],
        [1],
        [0],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0]])
        C=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [1,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0]
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        E=DriverInterval.bpcToE(b,q,C)
        t=20
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

    @staticmethod
    def motorTransmission1():
        A=Benchmarks.MotorTransmission1.A
        B=Benchmarks.MotorTransmission1.B
        mode='.'
        '''E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }'''
        b=np.array([
        [1],
        [0],
        [1],
        [0],
        [1],
        [1],
        [0]
        ])
        q=np.array([[mpi(-0.1,0.1),0,mpi(-0.1,0.1),0,0,mpi(-0.1,0.1),0]])
        C=np.array([
        [1,0,0,0,0,0,0],
        [0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0],
        [0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        ])
        E=DriverInterval.bpcToE(b,q,C)
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        U=[]
        E=DriverInterval.bpcToE(b,q,C)
        t=20
        vrfy=VerifyInterval(A,B,E,IS,t,U)
        print(vrfy.computeReachSet())
        print("")
        print(vrfy.computePerturbFreeReachSet())

class DriverRobustMetric:

    @staticmethod
    def illustExample():
        A=Benchmarks.IllustExample.A
        B=Benchmarks.IllustExample.B
        mode='.'
        t=20
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'

        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def stableSystem1():
        A=Benchmarks.StableSystem1.A
        B=Benchmarks.StableSystem1.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def stableSystem2():
        A=Benchmarks.StableSystem2.A
        B=Benchmarks.StableSystem2.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def stableSystem3():
        A=Benchmarks.StableSystem3.A
        B=Benchmarks.StableSystem3.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def stableSystem4():
        A=Benchmarks.StableSystem4.A
        B=Benchmarks.StableSystem4.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def flightEnvelope():
        A=Benchmarks.FlightEnvelope.A
        B=Benchmarks.FlightEnvelope.B
        mode='.'
        U=np.array([
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)],
        [(-5e6,-4e6)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def coOPVehiclesI():
        A=Benchmarks.CoOPVehiclesI.A
        B=Benchmarks.CoOPVehiclesI.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def coOPVehiclesII():
        A=Benchmarks.CoOPVehiclesII.A
        B=Benchmarks.CoOPVehiclesII.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def pkpd2():
        A=Benchmarks.PKPD2.A
        B=Benchmarks.PKPD2.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]


        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def spaceCraftRndzvs():
        A=Benchmarks.SpaceCraftRndzvs.A
        B=Benchmarks.SpaceCraftRndzvs.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def holesCXc():
        A=Benchmarks.HolesCXc.A
        B=Benchmarks.HolesCXc.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def motorTransmission1():
        A=Benchmarks.MotorTransmission1.A
        B=Benchmarks.MotorTransmission1.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())

    def motorTransmission2():
        A=Benchmarks.MotorTransmission2.A
        B=Benchmarks.MotorTransmission2.B
        mode='.'
        U=np.array([
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)],
        [(-2e5,-1e5)]
        ])
        IS=np.array([
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)],
        [mpi(1,3)]
        ])
        t=20
        method='Loan'
        rM=RobustMetric(A,B,U,method,IS,t,mode,0)
        print("Robustness Metric: ",rM.getMaxSafeBloat())







# Write your driver code Where
DriverRobustMetric.motorTransmission1()
