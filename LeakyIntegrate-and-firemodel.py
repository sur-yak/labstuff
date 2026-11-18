import numpy as np
import matplotlib.pyplot as plt

# Membrane parameters
Cm = 1       # nF
gL = 0.1       # µS
EL = -65.0     # mV

# Spike parameters
V_threshold = -50.0  # mV
V_reset = -65.0      # mV
t_ref = 2          # ms

# Simulation parameters
dt = 0.1       # ms
tmax = 1200    # ms
t = np.arange(0, tmax + dt, dt)

# Current injection
t0 = 100       # ms
t1 = 1100      # ms
Iext = np.zeros_like(t)

pulse = (t >= t0) & (t < t1)
Iext[pulse] = 2.0  # nA

# Initialize membrane potential
V = np.zeros_like(t)
V[0] = EL

# Store spike times
spike_times = []

# TODO: Implement Euler's method to update the membrane potential.

def dVdt(V_now, I_now):
    return (I_now - gL * (V_now - EL)) / Cm

def euler(V_now, I_now):
        return V_now + dt * dVdt(V_now, I_now)


# TODO: Detect threshold crossings, record spike times, and reset the voltage.


#Loop Euler's method until tmax is reached
refractory_period = 0

for i in range(len(t)-1):

    # TODO: Hold the voltage at V_reset during the refractory period.
    if t[i + 1] < refractory_period: 
        V[i + 1] = V_reset

    else:
        V[i + 1] = euler(V[i], Iext[i])

        if V[i + 1] >= V_threshold:
            spike_times.append(t[i+1])
            V[i + 1] = V_reset
            refractory_period = t[i+1] + t_ref

# TODO: Plot membrane potential and injected current against time.
# Plot membrane potential and injected current
fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# Plot membrane potential
ax[0].plot(t, V, label="Membrane Potential")
ax[0].axhline(V_threshold, linestyle="--", label="Threshold")

xlowlimit = 0
xuplimit = tmax
# Mark spike times

ax[0].scatter(spike_times, [V_threshold] * len(spike_times), color="red", marker="o", zorder=3)

ax[0].set_ylabel("Voltage (mV)")
ax[0].set_title("Leaky Integrate-and-Fire Neuron")
ax[0].set_xlim(xlowlimit, xuplimit)
ax[0].legend()
ax[0].grid()

# Plot injected current
ax[1].plot(t, Iext, color="orange", label="Injected Current")
ax[1].set_xlabel("Time (ms)")
ax[1].set_ylabel("Current (nA)")
ax[1].set_xlim(xlowlimit, xuplimit)
ax[1].legend()
ax[1].grid()

plt.tight_layout()
plt.show()

# TODO: Calculate the firing rate during current injection and plot f-I curve.

firing_rate = len(spike_times) / ((t1-t0) / 1000)
print("The firing rate at", Iext[pulse][0], "nA is", firing_rate, "Hz")

def firing_rate_graph(Ix):
    I_values = np.zeros_like(t)
    I_values[pulse] = Ix

    V_values = np.zeros_like(t)
    V_values[0] = EL

    spike_timesx = []
    refractory_periodx = 0

    for i in range(len(t)-1):
        if t[i+1] < refractory_periodx:
            V_values[i+1] = V_reset

        else:
            V_values[i+1] = euler(V_values[i], I_values[i])

            if V_values[i+1] >= V_threshold:
                spike_timesx.append(t[i + 1])
                V_values[i + 1] = V_reset
                refractory_periodx = t[i + 1] + t_ref
    return V_values, spike_timesx

currents = np.arange(0, 3.5, 0.5)
firing_rates = []
for i in currents:
    V_values, spike_timesx = firing_rate_graph(i)
    firing_ratex = len(spike_timesx) / ((t1-t0)/1000)
    firing_rates.append(firing_ratex)

plt.plot(currents, firing_rates)

plt.xlabel("Injected Current (nA)")
plt.ylabel("Firing Rate (Hz)")
plt.title("Firing Rate vs. Injected Current")

plt.grid()
plt.show()
