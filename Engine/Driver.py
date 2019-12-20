'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

from VerificationEngine import *
from Benchmarks import *

class Driver:

    @staticmethod
    def illustExample():
        A=Benchmarks.IllustExample.A
        B=Benchmarks.IllustExample.B
        mode='.'
        E={
        (0,0): [-2,2],
        (2,5): [-1,1],
        (4,6): [-2,2],
        (5,3): [-3,3]
        }
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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0.1,10,0.01)

    @staticmethod
    def flightEnvelope():
        A=Benchmarks.FlightEnvelope.A
        B=Benchmarks.FlightEnvelope.B
        mode='.'
        E={
        (3,7): [-1,1],
        (4,6): [-1,1],
        }
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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0,5,0.01)

    @staticmethod
    def pkpd2():
        A=Benchmarks.PKPD2.A
        B=Benchmarks.PKPD2.B
        mode='.'
        E={
        (0,4): [-1,1],
        }
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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0.1,3,0.01)

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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0.50,10


        ,0.01)

    @staticmethod
    def spaceCraftRndzvs():
        A=Benchmarks.SpaceCraftRndzvs.A
        B=Benchmarks.SpaceCraftRndzvs.B
        mode='.'
        E={
        (2,0): [-1,1],
        (2,3): [-0.1,0.1],
        (3,2): [-0.5,0.5],
        }
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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0,10,0.01)

    @staticmethod
    def holesCXc():
        A=Benchmarks.HolesCXc.A
        B=Benchmarks.HolesCXc.B
        mode='.'
        E={
        (0,3): [-0.5,0.5],
        (1,2): [-0.1,0.1],
        (3,2): [-0.005,0.005],
        (4,3): [-0.0045,0.0045]
        }
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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0,10,0.01)

    @staticmethod
    def holesCXc():
        A=Benchmarks.HolesCXc.A
        B=Benchmarks.HolesCXc.B
        mode='.'
        E={
        (0,3): [-0.5,0.5],
        (1,2): [-0.1,0.1],
        (3,2): [-0.005,0.005],
        (4,3): [-0.0045,0.0045]
        }
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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0,5,0.01)

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

        vrfy=Verify(A,B,E,IS,t,U)
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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0,3,0.01)

    @staticmethod
    def motorTransmission1():
        A=Benchmarks.MotorTransmission1.A
        B=Benchmarks.MotorTransmission1.B
        mode='.'
        E={
        (0,6): [-0.1,0.1],
        (1,6): [-1,1]
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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0,5,0.01)

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

        vrfy=Verify(A,B,E,IS,t,U)
        vrfy.plotTime(0,5,0.01)


# Write your driver code Where

Driver.flightEnvelope()
