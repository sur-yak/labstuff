# Libraries: matplotlib, numpy, pandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#setting up values
lam = 1
ome = 1
a = 1
b = 1
t = 0.01



def f(x, y, r):
    dxdt = (lam - b*(r**2))*x - (ome + a*(r**2))*y
    return dxdt

def g(x, y, r):
    dydt = (ome + a*(r**2))*x + (lam - b*(r**2))*y
    return dydt

#RK4
def RK4(x, y):
    #1st step
    r = x**2 + y**2
    k1 = f(x, y, r)
    l1 = g(x, y, r)
    x1 = x + (t/2)*k1
    y1 = y + (t/2)*l1
    r1 = (x1)**2 + (y1)**2

    #2nd step
    k2 = f(x1, y1, r1)
    l2 = g(x1, y1, r1)
    x2 = x + (t/2)*k2
    y2 = y + (t/2)*l2
    r2 = x2**2 + y2**2

    #3rd step
    k3 = f(x2, y2, r2)
    l3 = g(x2, y2, r2)
    x3 = x + t*k3
    y3 = y + t*l3
    r3 = (x3)**2 + (y3)**2

    #4th step
    k4 = f(x3, y3, r3)
    l4 = g(x3, y3, r3)
    x = x + (t/6)*(k1 + 2*(k2) + 2*(k3) + k4)
    y = y + (t/6)*(l1 + 2*(l2) + 2*(l3) + l4)
    return x, y

#maximum time
Tmax = 10 * np.pi

#setting up x-axis intervals and thier y-axis values
time_intervals = []
x_values = []
y_values = []
xx_values = []
yy_values = []

t_now = 0

#initial values for x and y at t=0
x = 1
y = 0
xx = np.sin(2*t+np.pi/2)
yy = np.cos(2*t+np.pi/2)

while t_now <= Tmax:

    #storing the values in list
    time_intervals.append(t_now)
    x_values.append(x)
    y_values.append(y)
    xx_values.append(xx)
    yy_values.append(yy)


    #assigning next RK4 step to variable
    x, y = RK4(x, y)
    xx, yy = RK4(xx, yy)

    #move to next interval
    t_now = t_now + t



#x vs t graph
plt.plot(time_intervals, x_values)
plt.plot(time_intervals, xx_values)
plt.xlabel("t")
plt.ylabel("x")
plt.title("x vs. t")

plt.show()

#y vs t graph
plt.plot(time_intervals, y_values)
plt.plot(time_intervals, yy_values)
plt.xlabel("t")
plt.ylabel("y")
plt.title("y vs. t")
    
plt.show()