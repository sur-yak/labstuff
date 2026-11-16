#Libraries: matplotlib, numpy, pandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#setting arrays and constants. Each position in the array corresponds to one oscillator

ome = np.array([1, 1])
a = np.array([1, 1])
lam = np.array([1, 2])
b = np.array([1, 2])

t = 0.01

holdvalue = 1.8

#r applies to x value
alpha1_r = np.array([
    [-1, -1],
    [-1, -1]
])

#imaginary value (haven't used it yet)
alpha1_i = np.array([
    [0, 0],
    [0, 0]
])

#alpha2 is written as beta in literature
alpha2_r = np.array([
    [-1, -1],
    [-1, -1]
])

alpha2_i = np.array([
    [0, 0],
    [0, 0]
])
tconn = 10.21
#RK4
def RK4(x, y, t_now):
    # setting up function to call for equation of each step
    def ddt(x, y, t_now):
        r2 = x**2 + y**2
        
        if t_now < tconn:
            A = np.zeros_like(x)
            B = np.zeros_like(y)

        else:
            A = alpha1_r @ x - alpha1_i @ y
            B = alpha2_i @ x + alpha2_r @ y

        dxdt = (lam - b*(r2))*x - (ome + a*(r2))*y + A
        dydt = (ome + a*(r2))*x + (lam - b*(r2))*y + B
        return dxdt, dydt

    #RK4 step A
    ka, la = ddt(x, y, t_now)

    #RK4 step B
    xa = x + (t/2)*ka
    ya = y + (t/2)*la
    kb, lb = ddt(xa, ya, t_now + t/2)

    #RK4 step C
    xb = x + (t/2)*kb
    yb = y + (t/2)*lb
    kc, lc = ddt(xb, yb, t_now + t/2)

    #RK4 step D
    xc = x + (t)*kc
    yc = y + (t)*lc
    kd, ld = ddt(xc, yc, t_now + t)

    x = x + (t/6)*(ka + 2*(kb) + 2*(kc) + kd)
    y = y + (t/6)*(la + 2*(lb) + 2*(lc) + ld)

    return x, y

Tmax = 75.5

#initial values for x and y at t=0
t_now = 0

x = np.array([1, 1])

y = np.array([0, 0])

#storage list, they store multiple arrays
time_intervals = []
x_values = []
y_values = []

x_initial = x.copy()
y_initial = y.copy()

currentclosestvaluetokeepat = 75.5

while t_now <= Tmax:
    time_intervals.append(t_now)

    x_values.append(x.copy())
    y_values.append(y.copy())

    if t_now < holdvalue:
        x, y = RK4(x, y, t_now)
        x[1] = x_initial[1]
        y[1] = y_initial[1]
    else:
        x, y = RK4(x, y, t_now)

    t_now = t_now + t

#convert x and y values to arrays
x_values = np.array(x_values)
y_values = np.array(y_values)

#plotting x vs. t
for i in range(len(x)):
    plt.plot(time_intervals, x_values[:, i], label=f"x{i+1}")

plt.xlabel("t")
plt.ylabel("Value")
plt.title("Oscillator x values vs. time")
plt.ylim(-2, 2)
plt.xlim(0, 75.5)
plt.axvline(x=tconn, color='red', linestyle=':', label=f'Tconn = {tconn}')
plt.legend()
plt.grid()
plt.show()