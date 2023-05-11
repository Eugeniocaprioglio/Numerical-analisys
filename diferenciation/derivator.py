import math
import numpy as np
import matplotlib.pyplot as plt

def fun(x):
    return x*math.exp(-x**2)

def derivator(f, x, dx):
    return ((f(x + dx) - f(x - dx))/(2*dx))
    
d = 0.0001
a = -5
b = 5
n = 1000

t = np.linspace(a, b, n)
y = np.zeros(n)
dy = np.zeros(n)
for i in range(n):
    y[i] = fun(t[i])
    dy[i] = derivator(fun, t[i], d)


fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.plot(t, y, color = color, label='function')
fig.set_label('time (s)')
ax1.set_ylabel('f(t)', color = color)

ax2 = ax1.twinx()
color = 'tab:green'
ax2.set_ylabel("f'(t)", color = color)

ax2.plot(t, dy, color = color, label='derivative')
fig.legend()
plt.grid(color=color)
plt.show()
