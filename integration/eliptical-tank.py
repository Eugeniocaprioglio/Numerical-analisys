import math
"""
    This program calculates the volume of an horizontal cylindrical-elliptical tank given the hight of the liquid.
    It also calculates the height of the liquid given the volume of liquid inside.
"""

# Tank meassurments

W = 2.3
H = 1.62
L = 7

# Total volume of the tank

max_vol = math.pi*W*H*L/4

d = .0001 # step size

# function used for calculate the area of the ellipse. Is valid for 0 < x < H

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
