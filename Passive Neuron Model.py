import numpy as np
import matplotlib.pyplot as plt

# Membrane parameters
#membrane capacitence
Cm = 1.0       # nF
#leak conducatnce
gL = 0.1       # µS
#leak reversal potential
EL = -65.0     # mV

# Simulation parameters
dt = 0.1       # ms
tmax = 500     # ms
t = np.arange(0, tmax + dt, dt)

# Current injection
t0 = 100       # ms
t1 = 300       # ms
Iext = np.zeros_like(t)

pulse = (t >= t0) & (t < t1)
Iext[pulse] = 4.0  # nA

# Initialize membrane potential
V = np.zeros_like(t)
V[0] = EL

# TODO: Implement Euler's method here
def euler(V_now, I_now):
    def dVdt(V_now, I_now):
        dVdt = (I_now - gL * (V_now - EL)) / Cm
        return dVdt
    
    V = V_now + dt * dVdt(V_now, I_now)

    return V

#Loop Euler's method until tmax is reached
for i in range(len(t)-1):
    V[i + 1] = euler(V[i], Iext[i])

# TODO: Plot membrane potential against time
plt.plot(t, V)
plt.ylim(-80, -40)
plt.xlim(0, tmax)
plt.xlabel("time (ms)")
plt.ylabel("Membrane Potentional (mV)")
plt.title("Passive Neuron Membrane Potential over Time")
plt.axhline(EL, color='purple', linestyle=':', label="Leaky Reversal Potential")
plt.axvline(t0, color='red', linestyle=':', label="current input on")
plt.axvline(t1, color='red', linestyle=':', label="current input off")
plt.legend()
plt.grid(True)
plt.show()