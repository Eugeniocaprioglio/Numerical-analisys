import math	
import numpy as np
import matplotlib.pyplot as plt



#------------Define the constants and calculate Fourier number---

#----dimensions of the discretization elements-----

delta_x = 0.01 # meters
delta_t = 0.15 # seconds
thickness = 0.01 # meters

#----material data------

molar_heat_capacity = 25418 # J/(kmol*K)
molar_mass = 197 # kg/kmol
mass_density = 19500 # kg/m^3
thermal_conductivity = 418 # W/(m*K)
mass_heat_capacity =  molar_heat_capacity/molar_mass # J/(kg*K)
alpha = thermal_conductivity/(mass_density*mass_heat_capacity) # m^2/s

def calculate_fourier():
    dt_max = 0.25*delta_x**2/alpha
    if (delta_t > dt_max):
        return 'The system is unstable. Choose a value of delta_T smaller than ' + str(round(dt_max, 4)) + ' sec or bigger value of delta_X'
    return alpha*delta_t/delta_x**2

#-----This function creates the initial matrix with different temperatures --
# ----in the inner squares----------------------------------------------------

def create_matrix(columns, raws): 

    mat = np.full((raws, columns),0)
    for i in range(int(columns/3),columns - int(columns/3)):
        for j in range(int(raws/3),raws - int(raws/3)):
            mat[j][i] = 100

    return np.array(mat)

def matrix_average(m):
        s = 0
        for i in range(1,len(m)-1):
            for j in range(1,len(m[0])-1):
                s += m[i,j]
        
        return s / ((len(m) - 1)*(len(m) - 1))

#-------Border conditions----------------------------------------------------------------------

def border_conditions(t):

    t[:,0] = 0 # t[:,1]  #left
    t[:,-1] = 0  # t[:,-2]    #right
    t[0,:] = 0 # t[1,:]   #bottom
    t[-1,:] = 200 #t[-2,:]     #top
    return t

 
def thermal_energy(x_length, y_length, iterations):
    mi = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    if isinstance(fourier, float):
        mf = heat_equation(mi, fourier, iterations)
    else:
        print(fourier)
        return
    c = mass_heat_capacity * mass_density * thickness * x_length * y_length * delta_x**2
    ei = c * matrix_average(mi)
    ef = c * matrix_average(mf)
    power = (ef - ei)/(delta_t * iterations)
    return 'Initial energy: ' + str(round(ei)) + 'J \n Final energy: ' + str(round(ef)) + 'J \n Power: ' + str(round(power))

# ------- This function calculates the laplacian part of the heat equation. -------------

def laplacian(t):
    
    up_v, down_v, right_v, left_v = [np.zeros((len(t), len(t[0]))) for i in range(4)]
    t = np.array(t)
    for i in range(len(t) - 1):
        up_v[i,:] = t[i+1,:]
        down_v[i+1,:] = t[i,:]
    for i in range(len(t[0]) - 1):
        right_v[:,i] = t[:,i+1]
        left_v[:,i+1] = t[:,i]
    return np.array(left_v) + np.array(right_v) + np.array(up_v) + np.array (down_v)

# ------ Solve the heat equation-------------------

def heat_equation(u, fo, k):
    for i in range(k):
        u = (1 - 4*fo)*np.array(u) + fo*np.array(laplacian(u))   
        u = border_conditions(u)
    return u

#------- This function calculates the temperature in a position vs time------

def temp_vs_time(u, fo, x_coord, y_coord, k):
    arr = []
    for i in range(k):
        u = border_conditions(u)
        u = (1 - 4*fo)*np.array(u) + fo*np.array(laplacian(u))
        arr.append([i*delta_t/60, u[x_coord, y_coord]])      
    return np.array(arr)

#------------Plotting functions---------------------------------

def plot_grid(x_length, y_length, iterations):
    mat = create_matrix(x_length, y_length)
    mat = border_conditions(mat)
    fourier = calculate_fourier()
    if isinstance(fourier, float):
        t = str(round(delta_t*iterations, 2))
        dx = str(round(x_length*delta_x, 2))
        dy = str(round(y_length*delta_x, 2))
        plt.pcolormesh(heat_equation(mat, fourier, iterations), cmap='inferno', vmax=200, vmin=0)
        plt.title('Temperature Distribution for time = ' + t + 's' + '\n' + 'height: ' + dy + 'm, width: ' + dx + 'm')
        plt.colorbar()
        plt.show()
    else:
        print(fourier)
    return

def plot_yprofile(x_length, y_length, iterations, x_coord):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    axis = plt.subplot()
    if isinstance(fourier, float):
        t = str(round(delta_t*iterations, 2))
        dx = str(round(x_length*delta_x, 2))
        dy = str(round(y_length*delta_x, 2))
        x = str(round(delta_x*x_coord, 2))
        solution = heat_equation(mat, fourier, iterations)
        axis.plot([i*delta_x for i in range(len(solution[:,x_coord]))], solution[:,x_coord])
        axis.set_xlabel('Y-coordinate')
        axis.set_ylabel('Temperature (Celsius)')
        plt.title('Y axis temperature profile for time = ' + t + ' s and x = ' + x + 'm. \n Width = ' + dx + ' m. Height = ' + dy +' m')
        plt.show()
    else:
        print(fourier)
    return

def plot_xprofile(x_length, y_length, iterations, y_coord):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    axis = plt.subplot()
    if isinstance(fourier, float):
        t = str(round(delta_t*iterations, 2))
        dx = str(round(x_length*delta_x, 2))
        dy = str(round(y_length*delta_x, 2))
        y = str(round(delta_x*y_coord, 2))
        solution = heat_equation(mat, fourier, iterations)
        axis.plot([i*delta_x for i in range(len(solution[y_coord,:]))],solution[y_coord,:])
        axis.set_xlabel('X-coordinate')
        axis.set_ylabel('Temperature (Celsius)')
        plt.title('X axis temperature profile for time = ' + t + ' s and y = ' + y + 'm. \n Width = ' + dx + ' m. Height = ' + dy +' m')
        plt.show()
    else:
        print(fourier)
    return

def plot_3d_surface(x_length, y_length, iterations):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    axis = plt.axes(projection='3d')
    if isinstance(fourier, float):
        solution = heat_equation(mat, fourier, iterations)
        x = np.arange(0, x_length, 1)
        y = np.arange(0, y_length, 1)
        Y, X = np.meshgrid(x, y)
        Z = solution[Y, X]
        axis.plot_surface(X*delta_x, Y*delta_x, Z, cmap='inferno')
        axis.set_xlabel('X')
        axis.set_ylabel('Y')
        axis.set_zlabel('Temperature (Celcius)')
        plt.show()
    else:
        print(fourier)
    return

def plot_contour(x_length, y_length, iterations):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    axis = plt.subplot()
    if isinstance(fourier, float):
        solution = heat_equation(mat, fourier, iterations)
        length = max(x_length, y_length)
        x = np.arange(0, length, 1)
        y = np.arange(0, length, 1)
        Y, X = np.meshgrid(x, y)
        Z = solution[Y,X]
        cs = axis.contour(X*delta_x, Y*delta_x, Z)
        axis.set_xlabel('X')
        axis.set_ylabel('Y')
        axis.clabel(cs, inline=True, fontsize=10)
        axis.set_title('Temperature contour')
        plt.show()
    else:
        print(fourier)
    return

def plot_temp_vs_time(x_length, y_length, x_coord, y_coord, iterations):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    axis = plt.subplot()
    if isinstance(fourier, float):
        solution = temp_vs_time(mat, fourier,  x_coord, y_coord, iterations)
        axis.plot(solution[:,0], solution[:,1])
        t = str(round(iterations*delta_t/60, 2))
        dx = str(round(x_length*delta_x, 2))
        dy = str(round(y_length*delta_x, 2))
        x = str(round(delta_x*x_coord, 2))
        y = str(round(delta_x*y_coord, 2))
        plt.title('Temperature vs time for x = ' + x + 'm and y = ' + y + ' m \n Final time : ' + t + ' minutes \n Width = ' + dx + ' m. Height = ' + dy +' m')
        axis.set_xlabel('Time (minutes)')
        axis.set_ylabel('Temperature (Celcius)')
        plt.show()
    else:
        print(fourier)
    return

plot_temp_vs_time(60,60,30,40,6000)