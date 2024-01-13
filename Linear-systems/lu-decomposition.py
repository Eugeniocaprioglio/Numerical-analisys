import math
import numpy as np 
from matplotlib import pyplot as plt 

#--this function solves a lower triangular system by forward substitution
def lower_triang(mti, b):

    x = np.zeros(len(b))
    for k in range(len(mti)):
        sum = 0
        for i in range(k):
            sum += mti[k][i]*x[i]
        x[k] = (b[k] - sum)/mti[k][k]       
    return x
def upper_triang(mts, b):
  
#--this function solves a upper triangular system by backward substitution
    x = np.zeros(len(b))
    
    for k in reversed(range(len(mts))):
        sum = 0
        for i in range(k,len(mts)):
            sum += mts[k][i]*x[i]
        x[k] = (b[k] - sum)/mts[k][k]
    return x

#--- this function returns the LU facorization matrices
def doolittle(a):
    print(a)
    n = len(a)
    l = np.identity(n)
    u = np.zeros((n,n))
    
    for k in range(n):
        for j in range(k,n):
            sum = 0
            for p in range(k):
                sum += l[k][p]*u[p][j]
            u[k][j] = a[k][j] - sum
        for i in range(k+1,n):
            sum = 0
            for p in range(k):
                sum += l[i][p]*u[p][k]
            l[i][k] = 1/u[k][k]*(a[i][k] - sum)
    return [l,u]
 # This function performs the gauss elimination metod in order to solve the system m*x = b
def gauss_elimination(m,b):
    l,u = doolittle(m)
    y = lower_triang(l,b)
    x = upper_triang(u,y)
    return x
  
