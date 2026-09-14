# Libraries: matplotlib, numpy, pandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def RK4(x, y, ome, phi, a, b, h):
    # function definitions
    r = (phi/b)**0.5
    def f(x, y):
        return (phi - b*(r**2))*x - (ome + a*(r**2))*y
    def g(x, y):
        return (ome + a*(r**2))*x + (phi - b*(r**2))*y

    # 1st RK step
    k1 = -2*f(x, y)
    l1 = 2*g(x, y)

    # 2nd RK step
    k2 = -2*f(x + 0.5*h*l1, y + 0.5*h*k1)
    l2 = 2*g(x + 0.5*h*l1, y + 0.5*h*k1)

    # 3rd RK step
    k3 = -2*f(x + 0.5*h*l2, y + 0.5*h*k2)
    l3 = 2*g(x + 0.5*h*l2, y + 0.5*h*k2)

    # 4th RK step
    k4 = -2*f(x + h*l3, y + h*k3)
    l4 = 2*g(x + h*l3, y + h*k3)

    # final answer
    xnext = x + (h/6)*(l1 + 2*l2 + 2*l3 + l4)
    ynext = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    return xnext, ynext

# graphing RK4
x = 7
y = 4

x_values = [x]
y_values = [y]

for i in range(100):
    x, y = RK4(x, y, 1, 1, 1, 1, 0.1)
    x_values.append(x)
    y_values.append(y)

plt.plot(x_values, y_values)
plt.xlabel("x")
plt.ylabel("y")
plt.title("RK4 Method Test")
plt.grid()
plt.show()