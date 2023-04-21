import math


W = 2.3
H = 1.62
L = 7
max_vol = math.pi*W*H*L/4
d = .0001

def oval(x):
    return W*math.sqrt(1 - (2*x/H - 1)**2)


def trapezoid(f, a, b, max_area, dx):
    k = int((b-a)/dx)
    s = 0
    for i in range(k-1):
        s += f(a + (i+1)*dx)
        if (f(a) + 2*s)*dx/2 > max_area:
            print("The value of x for the max area is " + str(round(a + i*dx,4)))
            return  a + i*dx
        
    return (f(a) + 2*s + f(b))*dx/2

def volume(h):
    return trapezoid(oval, 0, h, 2*max_vol/L, d)*L

def liquid_height(v):
    if v > max_vol:
        return "The maximum volume of the tank is " + str(round(max_vol, 3))
    return trapezoid(oval, 0, H, v/L, d)
