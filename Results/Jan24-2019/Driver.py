'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''
from Benchmarks import *
from RevEigenAPI import *
from ReachSetAPI import *
from mpmath import *

class DriverDec:

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
    def stableSystem1():
        A=Benchmarks.StableSystem1.A
        B=Benchmarks.StableSystem1.B
        mode='.'
        dVal={0:[-0.01,0.01]}
        dVec={
        (0,0): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        ])
        t=5
        print("Stable System 1")
        ReachSet(A,B,IS,dVec,dVal,t).printReport()
        print("\n********************\n")
        RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()

    def stableSystem2():
        A=Benchmarks.StableSystem2.A
        B=Benchmarks.StableSystem2.B
        dVal={0:[-0.01,0.01]}
        dVec={
        (1,0): [-0.01,0.01]
        }
        IS=np.array([
        [1],
        [1],
        ])
        t=5
        print("Stable System 2")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()
        print("\n********************\n")
        RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()

    def stableSystem3():
        A=Benchmarks.StableSystem3.A
        B=Benchmarks.StableSystem3.B
        mode='.'
        dVal={0:[-0.01,0.01], 1:[-0.02,0.02]}
        dVec={
        (0,1): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
        IS=np.array([
        [1],
        [1],
        ])
        t=5
        print("Stable System 3")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        ReachSet(A,B,IS,dVec,dVal,t).printReport()
        print("\n********************\n")
        RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()

    def stableSystem4():
        A=Benchmarks.StableSystem4.A
        B=Benchmarks.StableSystem4.B
        mode='.'
        dVal={0:[-0.01,0.01], 1:[-0.02,0.02]}
        dVec={
        (0,1): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
        IS=np.array([
        [1],
        [1],
        [1]
        ])
        t=5
        print("Stable System 4")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        ReachSet(A,B,IS,dVec,dVal,t).printReport()
        print("\n********************\n")
        RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()

    @staticmethod
    def illustExample():
        A=Benchmarks.IllustExample.A
        B=Benchmarks.IllustExample.B

    @staticmethod
    def flightEnvelope():
        A=Benchmarks.FlightEnvelope.A
        B=Benchmarks.FlightEnvelope.B
        mode='.'
        dVal={0:[-0.01,0.01], 1:[-0.02,0.02]}
        dVec={
        (0,1): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
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
        t=5
        print("Flight Envelope")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()

    @staticmethod
    def coOPVehiclesI():
        A=Benchmarks.CoOPVehiclesI.A
        B=Benchmarks.CoOPVehiclesI.B
        mode='.'
        dVal={0:[-0.01,0.01], 1:[-0.02,0.02]}
        dVec={
        (0,1): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
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
        t=5
        print("CoOP Vehicles I")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()

    def coOPVehiclesII():
        A=Benchmarks.CoOPVehiclesII.A
        B=Benchmarks.CoOPVehiclesII.B
        dVal={0:[-0.01,0.01], 1:[-0.02,0.02]}
        dVec={
        (0,1): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
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
        t=5
        print("CoOP Vehicles II")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()

    @staticmethod
    def pkpd2():
        A=Benchmarks.PKPD2.A
        B=Benchmarks.PKPD2.B
        mode='.'
        dVal={0: [-0.01,0.01]}
        dVec={
        (0,1): [-0.1,0.1]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        t=5
        print("PKPD2")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()
        print("\n********************\n")
        RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()

    @staticmethod
    def dcConv():
        A=Benchmarks.DCConv.A
        B=Benchmarks.DCConv.B
        mode='.'
        dVal={0: [-0.01,0.01]}
        dVec={
        (0,1): [-0.1,0.1]
        }
        IS=np.array([
        [1],
        [1],
        [1]
        ])
        t=5
        print("DC Converter")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        ReachSet(A,B,IS,dVec,dVal,t).printReport()
        print("\n********************\n")
        RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()

    @staticmethod
    def spaceCraftRndzvs():
        A=Benchmarks.SpaceCraftRndzvs.A
        B=Benchmarks.SpaceCraftRndzvs.B
        mode='.'
        dVal={}
        dVec={}
        print("Space Craft Rndzvs")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        dVal={0:[-0.01,0.01], 1:[-0.02,0.02]}
        dVec={
        (0,1): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        t=5
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()

    @staticmethod
    def holesCXc():
        A=Benchmarks.HolesCXc.A
        B=Benchmarks.HolesCXc.B
        mode='.'
        dVal={0:[-0.01,0.01], 1:[-0.02,0.02]}
        dVec={
        (0,1): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
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
        t=5
        print("Holes")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()

    @staticmethod
    def holesPDp():
        A=Benchmarks.HolesPDp.A
        B=Benchmarks.HolesPDp.B

    @staticmethod
    def holesPXp():
        A=Benchmarks.HolesPXp.A
        B=Benchmarks.HolesPXp.B

    @staticmethod
    def motorTransmission1():
        A=Benchmarks.MotorTransmission1.A
        B=Benchmarks.MotorTransmission1.B
        dVal={0:[-0.01,0.01], 1:[-0.02,0.02]}
        dVec={
        (0,1): [-0.1,0.1],
        (1,0): [-0.01,0.01]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        t=5
        print("Motor Transmission 1")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()

    @staticmethod
    def motorTransmission2():
        A=Benchmarks.MotorTransmission2.A
        B=Benchmarks.MotorTransmission2.B
        mode='.'
        dVal={0:[-0.01,0.01]}
        dVec={
        (0,1): [-0.01,0.01],
        (1,0): [-0.01,0.01]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        t=5
        print("Motor Transmission 2")
        #r=RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()
        r=ReachSet(A,B,IS,dVec,dVal,t).printReport()
        print("\n********************\n")
        RevEigenDecomp(DriverDec.createMatrix(A,B,'.',0),dVal,dVec).printReport()







# Write your driver code Where
DriverDec.pkpd2()

