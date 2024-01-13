import math
import numpy as np 

def errore(t):
    sum = 0
    for i in range(len(t)):
        sum += t[i]**2
    return math.sqrt(sum)

def step_jacobi(matrix, x_0, b):
    n = len(x_0)
    x = np.zeros(n)
    for i in range(n):
        sum = 0
        for j in range(n):
            if j != i:
                sum += matrix[i][j]*x_0[j]
        x[i] = 1/matrix[i][i]*(b[i] - sum)
    return np.array(x)

def met_jacobi(a,b,x0,n_max, toll):
    sol1 = np.array(x0)
    sol2 = np.array(x0)
    d = [1 for i in range(len(sol1))]
    for i in range(n_max):
        sol2 = step_jacobi(a,sol1,b)
        d = sol2 - sol1
        sol1 = sol2
        if errore(d) < toll:
            break 
    return sol1
