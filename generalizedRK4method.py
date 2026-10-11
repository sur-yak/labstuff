#Libraries: matplotlib, numpy, pandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#setting arrays and constants. Each position in the array corresponds to one oscillator

lam = np.array([1, 1, 1])
ome = np.array([1, 1, 1])
a = np.array([1, 1, 1])
b = np.array([1, 1, 1])

t = 0.01

#r applies to x value
alpha_r = np.array([
    [0, -1, 0],
    [-1, 0, 0],
    [0, 0, -1]
])

#i applies to y value to make it 2-dimensional
alpha_i = np.array([
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
])

#RK4
def RK4(x, y):
    # setting up function to call for equation of each step
    def ddt(x, y):
        r = x**2 + y**2
        A = alpha_r @ x - alpha_i @ y
        B = alpha_i @ x + alpha_r @ y

        dxdt = (lam - b*(r**2))*x - (ome + a*(r**2))*y + A
        dydt = (ome + a*(r**2))*x + (lam - b*(r**2))*y + B
        return dxdt, dydt

    #RK4 step A
    ka, la = ddt(x, y)

    #RK4 step B
    xa = x + (t/2)*ka
    ya = y + (t/2)*la
    kb, lb = ddt(xa, ya)

    #RK4 step C
    xb = x + (t/2)*kb
    yb = y + (t/2)*lb
    kc, lc = ddt(xb, yb)

    #RK4 step D
    xc = x + (t)*kc
    yc = y + (t)*lc
    kd, ld = ddt(xc, yc)

    x = x + (t/6)*(ka + 2*(kb) + 2*(kc) + kd)
    y = y + (t/6)*(la + 2*(lb) + 2*(lc) + ld)

    return x, y

Tmax = 10*np.pi

#initial values for x and y at t=0
t_now = 0

x = np.array(
    [1, 2, 1.5]
)

y = np.array(
    [0.9, 0.4, 0.2]
)

#storage list, they store multiple arrays
time_intervals = []
x_values = []
y_values = []

while t_now <= Tmax:
    time_intervals.append(t_now)

    x_values.append(x.copy())
    y_values.append(y.copy())

    x, y = RK4(x, y)

    t_now = t_now + t

#convert x and y values to arrays
x_values = np.array(x_values)
y_values = np.array(y_values)

#plotting x vs t
for i in range(len(x)):
    plt.plot(time_intervals, x_values[:, i], label=f"Oscillator {i+1}")

plt.xlabel("t")
plt.ylabel("x")
plt.title("x vs. t")
plt.legend()
plt.grid()

plt.show()
#plotting y vs t
for i in range(len(y)):
    plt.plot(time_intervals, y_values[:, i], label=f"Oscillator {i+1}")

plt.xlabel("t")
plt.ylabel("y")
plt.title("y vs. t")
plt.legend()
plt.grid()

plt.show()

