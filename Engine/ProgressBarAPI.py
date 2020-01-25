from time import sleep
import sys

class ProgressBar:
    def __init__(self,max,min,step)

    def showProgress(self):

        for i in range(self.min,self.max):
            if (i%self.)
            sys.stdout.write('\r')
            # the exact output you're looking for:
            sys.stdout.write("[%-20s] %d%%" % ('='*i, 5*i))
            sys.stdout.flush()
            sleep(0.25)
