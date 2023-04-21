import numpy as np
import math
import matplotlib.pyplot as plt

#---------------Runge Kutta single step function-------------

def rk4_step(f, dt, tn, yn):
    k1 = f(tn, yn)
    k2 = f(tn + dt/2, yn + k1*dt/2)
    k3 = f(tn + dt/2, yn + 0.5*k2*dt)
    k4 = f(tn + dt, yn + k3*dt)
    return  yn + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)



#-----Define the differential equation-----------------------

def projectile(t, s):

    x, y, z, vx, vy, vz = s
    wx, wy = w
    return np.array([
        vx,
        vy,
        vz,
        -(mu/m)*(vx - wx)*pow((vx - wx)**2 + (vy - wy)**2 + vz**2, (p - 1)/2),
        -(mu/m)*(vy - wy)*pow((vx - wx)**2 + (vy - wy)**2 + vz**2, (p - 1)/2),
        -g - (mu/m)*(vz)*pow((vx - wx)**2 + (vy - wy)**2 + vz**2 , (p - 1)/2)
    ])


#  Initial position
x0 = 0
y0 = 0
z0 = 0

#  Initial speed and shooting angles
v0 = 80
theta = 15 #---azimut
phi = 30  #----inclination

#  Vector coordinates of the initial velocity
v0x = v0*math.cos(math.radians(theta))*math.cos(math.radians(phi))
v0y = v0*math.sin(math.radians(theta))*math.cos(math.radians(phi))
v0z = v0*math.sin(math.radians(phi))

#   Initial conditions array
initial_conditions = [x0, y0, z0, v0x, v0y, v0z]

#------Set the constants of the model-------------------------

g = 9.81 #gravitational intensity near the surfice of the rocky celestial body

#----Wind speed and direction-------
wind_speed = 0  #---m/s
wind_direction = 70 #--degrees with respect of the positive x-semiaxis (plane y = 0)
w = [wind_speed*math.cos(math.radians(wind_direction)), 
        wind_speed*math.sin(math.radians(wind_direction))] #--[wx, wy]


#---mass of the ball
m = 0.046

#---diameter of the ball
diameter = 0.042

#----drag force parameters------------
cd = 0.42 #-drag coefficient for a sphere
fluid_density = 1.25 #kg/m^3
area = 3.1416 * (diameter/2)**2
mu = 0.5 * cd * area * fluid_density
p = 1.8


def impact_floor(v, floor):
    if v[2] < floor:
        return True
    return False

#Compute trajectory

dt = 0.01
T = 15
time_points = int(T/dt)
t = np.linspace(0, T, time_points)

Y = np.zeros((6, time_points))
Y[:,0] = initial_conditions
yin = initial_conditions
for i in range(time_points - 1):
    yout = rk4_step(projectile, dt, t[i], yin)
    Y[:, i+1] = yout
    yin = yout
    floor_pos = -0.001
    if impact_floor(yout, floor_pos):
        yin[5] = -0.6*yout[5]
    


#Plot trajectory

ax = plt.figure().add_subplot(projection='3d')
ax.plot(Y[0,:], Y[1,:], Y[2,:], 'b')
plt.show()
