#Libraries: matplotlib, numpy, pandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#setting up values
lam1 = 1
ome1 = 1
a1 = 1
b1 = 1

lam2 = 1
ome2 = 1
a2 = 1
b2 = 1

t = 0.01

#real values
alpha11r = -1
alpha12r = -1
alpha21r = -1
alpha22r = -1

#imaginary values
alpha11i = 0
alpha12i = 0
alpha21i = 0
alpha22i = 0

def f(x, y, r, A):
    dx1dt = (lam1 - b1*(r**2))*x - (ome1 + a1*(r**2))*y + A
    return dx1dt

def g(x, y, r, B):
    dy1dt = (ome1 + a1*(r**2))*x + (lam1 - b1*(r**2))*y + B
    return dy1dt

def h(x, y, r, A):
    dx2dt = (lam2 - b2*(r**2))*x - (ome2 + a2*(r**2))*y + A
    return dx2dt

def i(x, y, r, B):
    dy2dt = (ome2 + a2*(r**2))*x - (ome2 + a1*(r**2))*y + B
    return dy2dt

#RK41
def RK4(xa, ya, xb, yb):
    #1st step
    ra = xa**2 + ya**2
    rb = xb**2 + yb**2
    A1 = (alpha11r)*xa + (alpha12r)*xb + (alpha11i)*ya + (alpha12i)*yb
    A2 = (alpha21r)*xa + (alpha22r)*xb + (alpha21i)*ya + (alpha22i)*yb
    B1 = (alpha11i)*xa + (alpha12i)*xb + (alpha11r)*ya + (alpha12r)*yb
    B2 = (alpha21i)*xa + (alpha22i)*xb + (alpha21r)*ya + (alpha22r)*yb

    k1 = f(xa, ya, ra, A1)
    l1 = g(xa, ya, ra, B1)
    m1 = h(xb, yb, rb, A2)
    n1 = i(xb, yb, rb, B2)
    xa1 = xa + (t/2)*k1
    ya1 = ya + (t/2)*l1
    xb1 = xb + (t/2)*m1
    yb1 = yb + (t/2)*n1
    ra1 = (xa1)**2 + (ya1)**2
    rb1 = (xb1)**2 + (yb1)**2
    A11 = (alpha11r)*xa1 + (alpha12r)*xb1 + (alpha11i)*ya1 + (alpha12i)*yb1
    A21 = (alpha21r)*xa1 + (alpha22r)*xb1 + (alpha21i)*ya1 + (alpha22i)*yb1
    B11 = (alpha11i)*xa1 + (alpha12i)*xb1 + (alpha11r)*ya1 + (alpha12r)*yb1
    B21 = (alpha21i)*xa1 + (alpha22i)*xb1 + (alpha21r)*ya1 + (alpha22r)*yb1

    #2nd step
    k2 = f(xa1, ya1, ra1, A11)
    l2 = g(xa1, ya1, ra1, B11)
    m2 = h(xb1, yb1, rb1, A21)
    n2 = i(xb1, yb1, rb1, B21)
    xa2 = xa1 + (t/2)*k2
    ya2 = ya1 + (t/2)*l2
    xb2 = xb1 + (t/2)*m2
    yb2 = yb1 + (t/2)*n2
    ra2 = xa2**2 + ya2**2
    rb2 = xb2**2 + yb2**2
    A12 = (alpha11r)*xa2 + (alpha12r)*xb2 + (alpha11i)*ya2 + (alpha12i)*yb2
    A22 = (alpha21r)*xa2 + (alpha22r)*xb2 + (alpha21i)*ya2 + (alpha22i)*yb2
    B12 = (alpha11i)*xa2 + (alpha12i)*xb2 + (alpha11r)*ya2 + (alpha12r)*yb2
    B22 = (alpha21i)*xa2 + (alpha22i)*xb2 + (alpha21r)*ya2 + (alpha22r)*yb2

    #3rd step
    k3 = f(xa2, ya2, ra2, A12)
    l3 = g(xa2, ya2, ra2, B12)
    m3 = h(xb2, yb2, rb2, A22)
    n3 = i(xb2, yb2, rb2, B22)
    xa3 = xa2 + t*k3
    ya3 = ya2 + t*l3
    xb3 = xb2 + t*m3
    yb3 = yb2 + t*n3
    ra3 = (xa3)**2 + (ya3)**2
    rb3 = (xb3)**2 + (yb3)**2
    A13 = (alpha11r)*xa3 + (alpha12r)*xb3 + (alpha11i)*ya3 + (alpha12i)*yb3
    A23 = (alpha21r)*xa3 + (alpha22r)*xb3 + (alpha21i)*ya3 + (alpha22i)*yb3
    B13 = (alpha11i)*xa3 + (alpha12i)*xb3 + (alpha11r)*ya3 + (alpha12r)*yb3
    B23 = (alpha21i)*xa3 + (alpha22i)*xb3 + (alpha21r)*ya3 + (alpha22r)*yb3

    #4th step
    k4 = f(xa3, ya3, ra3, A13)
    l4 = g(xa3, ya3, ra3, B13)
    m4 = h(xb3, ya3, rb3, A23)
    n4 = i(xb3, ya3, rb3, B23)

    xa = xa + (t/6)*(k1 + 2*(k2) + 2*(k3) + k4)
    ya = ya + (t/6)*(l1 + 2*(l2) + 2*(l3) + l4)
    xb = xb + (t/6)*(m1 + 2*(m2) + 2*(m3) + m4)
    yb = yb + (t/6)*(n1 + 2*(n2) + 2*(n3) + n4)
    return xa, ya, xb, yb

#maximum time
Tmax = 10 * np.pi

#setting up x-axis intervals and thier y-axis values
time_intervals = []
xa_values = []
ya_values = []
xb_values = []
yb_values = []
xx_values = []
yy_values = []

t_now = 0

#initial values for x and y at t=0
xa = 1
ya = 0
xb = 1
yb = 0
xx = np.sin(2*t+np.pi/2)
yy = np.cos(2*t+np.pi/2)

while t_now <= Tmax:

    #storing the values in list
    time_intervals.append(t_now)
    xa_values.append(xa)
    ya_values.append(ya)
    xb_values.append(xb)
    yb_values.append(yb)
    xx_values.append(xx)
    yy_values.append(yy)


    #assigning next RK4 step to variable
    xa, ya, xb, yb = RK4(xa, ya, xb,yb)

    #move to next interval
    t_now = t_now + t



#x vs t graph
plt.plot(time_intervals, xa_values)
plt.plot(time_intervals, xb_values)
plt.xlabel("t")
plt.ylabel("x")
plt.title("x vs. t")

plt.show()

#y vs t graph
plt.plot(time_intervals, ya_values)
plt.plot(time_intervals, yb_values)
plt.xlabel("t")
plt.ylabel("y")
plt.title("y vs. t")
    
plt.show()