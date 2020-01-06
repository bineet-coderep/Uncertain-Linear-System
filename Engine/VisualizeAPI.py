'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

Given a vector of interval and the dimensions to visualize,
it plots a corresponding rectangle (to be extended to other shapes)
'''

import numpy as np
import matplotlib.pyplot as plt

class Visualize:
    '''
    Given a vector of interval and the dimensions to visualize,
    it plots a corresponding rectangle (to be extended to other shapes)
    '''

    def __init__(self,v,x,y,shape='rect'):
        self.vectr=v # the vector to visualize (polytope)
        self.x=x # x axis, given as an index of v
        self.y=y # y axis, given as an index of v
        self.n=v.shape[0] # dimension of the polytope
        self.shape=shape # only rectangle is allowed as of now

    def draw(self):
        '''
        Draws the self rectangle
        '''
        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        x_interval=self.vectr[self.x]
        y_interval=self.vectr[self.y]
        xy=(x_interval[0],y_interval[0])
        width=x_interval[1]-x_interval[0]
        height=y_interval[1]-y_interval[0]
        plt.gca().add_patch(plt.Rectangle(xy,width,height,0.0,fill=True,fc='#5832a8'))
        plt.axis('scaled')
        plt.show()

        '''plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)

        circle = plt.Rectangle((100,1000),10,10,0.0,fill=True,fc='#03bafc')
        circle2 = plt.Rectangle((2,2),5,5,0.0,fill=True,fc='#7303fc')

        plt.gca().add_patch(circle)
        #plt.gca().add_patch(circle2)

        plt.axis('scaled')
        plt.show()'''

    def drawCompare(self,s):
        '''
        Draws both self and given parameter s of type Visualize
        '''
        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        x_interval=self.vectr[self.x]
        y_interval=self.vectr[self.y]
        xy=(x_interval[0],y_interval[0])
        width=x_interval[1]-x_interval[0]
        height=y_interval[1]-y_interval[0]
        x2_interval=s.vectr[s.x]
        y2_interval=s.vectr[s.y]
        xy2=(x2_interval[0],y2_interval[0])
        width2=x2_interval[1]-x2_interval[0]
        height2=y2_interval[1]-y2_interval[0]
        if width2>width and height2>height:
            plt.gca().add_patch(plt.Rectangle(xy2,width2,height2,0.0,fill=True,fc='#7303fc'))
            plt.gca().add_patch(plt.Rectangle(xy,width,height,0.0,fill=True,fc='#03bafc'))
        else:
            plt.gca().add_patch(plt.Rectangle(xy,width,height,0.0,fill=True,fc='#03bafc'))
            plt.gca().add_patch(plt.Rectangle(xy2,width2,height2,0.0,fill=True,fc='#7303fc'))
        plt.axis('scaled')
        if self.x==s.x and self.y==s.y:
            plt.xlabel("Dimension "+str(self.x))
            plt.ylabel("Dimension "+str(self.y))
        else:
            plt.xlabel("Dimension Mismatch")
            plt.ylabel("Dimension Mismatch")
        plt.show()

if False:
    v=np.array([
    [1,1.2],
    [1,2],
    [1,2],
    [1,2],
    ])
    x=0
    y=2
    vs=Visualize(v,x,y)

    v2=np.array([
    [-1,3],
    [-1,3],
    [-1,3],
    [-1,3],
    ])
    x2=0
    y2=2
    vs2=Visualize(v2,x2,y2)
    vs.drawCompare(vs2)
