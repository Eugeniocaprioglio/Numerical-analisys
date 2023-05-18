
import numpy as np

d = 2
L = 3

def volume(h):
    return ((d/2)**2 * np.arccos((d/2-h)/(d/2)) - (d/2-h)*np.sqrt((d/2)**2 - (d/2 - h)**2))*L - 5

def newton(f, x0, dx, n):
    x = x0
    for i in range(n):
        x = x - f(x)/((f(x + dx) - f(x - dx))/(2*dx))

    return x



sol = newton(volume, 0.2, 0.0001, 5)
print(sol)
print(volume(sol))

