'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

from ComputeU import *
from Benchmarks import *

mode='.'
h=0.01

class DriverCompU:

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

    def stableSystem1():
        dynA=Benchmarks.StableSystem1.A
        dynB=Benchmarks.StableSystem1.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        E={
        (2,2): [0.95,1.05]
        }
        C=[0,0,0]
        V=np.array([
        [1,0,0],
        [0,1,0],
        [0,0,1],
        ])
        P=[(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def stableSystem2():
        dynA=Benchmarks.StableSystem2.A
        dynB=Benchmarks.StableSystem2.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E={
        (0,1): [1-(P/100),1+(P/100)]
        }
        C=[0,0]
        V=np.array([
        [1,0],
        [0,1]
        ])
        P=[(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def stableSystem3():
        dynA=Benchmarks.StableSystem3.A
        dynB=Benchmarks.StableSystem3.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E={
        (0,1): [1-(P/100),1+(P/100)]
        }
        C=[0,0]
        V=np.array([
        [1,0],
        [0,1]
        ])
        P=[(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def stableSystem4():
        dynA=Benchmarks.StableSystem4.A
        dynB=Benchmarks.StableSystem4.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E={
        (0,0): [1-(P/100),1+(P/100)],
        (2,2): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0]
        V=np.array([
        [1,0,0],
        [0,1,0],
        [0,0,1],
        ])
        P=[(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def flightEnvelope():
        dynA=Benchmarks.FlightEnvelope.A
        dynB=Benchmarks.FlightEnvelope.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E_weird={
        (3,7): [1-(P/100),1+(P/100)],
        (4,6): [1-(P/100),1+(P/100)]
        }
        E={
        (0,1): [1-(P/100),1+(P/100)],
        (0,3): [1-(P/100),1+(P/100)],
        (1,4): [1-(P/100),1+(P/100)],
        (8,11): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
        ])
        P=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def coOPVehiclesI():
        dynA=Benchmarks.CoOPVehiclesI.A
        dynB=Benchmarks.CoOPVehiclesI.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E={
        (0,0): [1-(P/100),1+(P/100)],
        (2,0): [1-(P/100),1+(P/100)],
        (2,1): [1-(P/100),1+(P/100)],
        (2,8): [1-(P/100),1+(P/100)],
        (5,0): [1-(P/100),1+(P/100)],
        (5,8): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,1]
        ])
        P=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def pkpd2():
        dynA=Benchmarks.PKPD2.A
        dynB=Benchmarks.PKPD2.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=20
        E={
        (0,4): [1-(P/100),1+(P/100)],
        (3,3): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        P=[(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def dcConv():
        dynA=Benchmarks.DCConv.A
        dynB=Benchmarks.DCConv.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E={
        (0,2): [1-(P/100),1+(P/100)],
        (1,1): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0]
        V=np.array([
        [1,0,0],
        [0,1,0],
        [0,0,1]
        ])
        P=[(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def spaceCraftRndzvs():
        dynA=Benchmarks.SpaceCraftRndzvs.A
        dynB=Benchmarks.SpaceCraftRndzvs.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (2,1): [1-(P/100),1+(P/100)],
        (2,0): [1-(P/100),1+(P/100)],
        (2,3): [1-(P/100),1+(P/100)],
        (3,2): [1-(P/100),1+(P/100)],
        }
        C=[0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0],
        [0,1,0,0,0,0],
        [0,0,1,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,1,0],
        [0,0,0,0,0,1],
        ])
        P=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def holesCXc():
        dynA=Benchmarks.HolesCXc.A
        dynB=Benchmarks.HolesCXc.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (0,3): [0.9,1.1],
        (1,2): [0.9,1.1],
        (3,2): [0.9,1.1],
        (4,3): [0.9,1.1]
        }
        C=[0,0,0,0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,1],
        ])
        P=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def motorTransmission1():
        dynA=Benchmarks.MotorTransmission1.A
        dynB=Benchmarks.MotorTransmission1.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E={
        (0,6): [1-(P/100),1+(P/100)],
        (1,6): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0],
        [0,1,0,0,0,0,0],
        [0,0,1,0,0,0,0],
        [0,0,0,1,0,0,0],
        [0,0,0,0,1,0,0],
        [0,0,0,0,0,1,0],
        [0,0,0,0,0,0,1]
        ])
        P=[(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)

    def motorTransmission2():
        dynA=Benchmarks.MotorTransmission2.A
        dynB=Benchmarks.MotorTransmission2.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (4,0): [0.9,1.1],
        (4,1): [0.9,1.1]
        }
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1]
        ])
        P=[(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        u=CompU(A,E)
        u.computeU(rs)



DriverCompU.motorTransmission2()
