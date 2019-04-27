import math
from decimal import *
import numpy as np

getcontext().prec = 100

#def f(x):
#    y = (Decimal(x[0][0])**2)+(Decimal(x[1][0])**2)
#    return y

def f(x):
    y = (x[0][0]-x[1][0])**4 + (x[0][0]+2*x[1][0]-3)**2
    return y

def gr(x):
    d1 = (4*(x[0][0]-x[1][0])**3)+(2*(x[0][0]+2*x[1][0]-3))
    d2 = (-4*(x[0][0]-x[1][0])**3)+(4*(x[0][0]+2*x[1][0]-3))
    gra = np.array([[d1],[d2]])
    return gra

def esi(x):
    d1 = (12*(x[0][0]-x[1][0])**2)+2
    d2 = (-12*(x[0][0]-x[1][0])**2)+4
    d3 = (-12*(x[0][0]-x[1][0])**2)+4
    d4 = (12*(x[0][0]-x[1][0])**2)+8
    esia = np.array([[d1,d2],[d3,d4]])
    return esia


def metodoNewton(x):
    while(True):
        g = gr(x)
        if(np.linalg.norm(g)<0.001):
            break
        e = esi(x)
        d = -(np.linalg.solve(e,g))
        #print d
        #print x
        x = x+d
        #print x
    return x



x = np.array([[Decimal(5)],[Decimal(10)]],dtype=np.dtype(Decimal))
#print f(x)
print (gr(x))
print ('---------')
print (esi(x))
print ('======>')
print (np.linalg.solve(esi(x),gr(x)))
#print metodoNewton(x)
