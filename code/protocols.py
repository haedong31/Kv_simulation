# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 2020

@author: Haedong Kim

Andrew's protocols
"""


import numpy as np
import matplotlib.pyplot as plt


## functions for generating protocols -----
def voltage_clamp(t, holding_p, holding_t, P1, P1_t, P2):
    v = np.piecewise(t, 
                     [t < holding_t, (t >= holding_t) & (t <= P1_t), t > P1_t],
                     [holding_p, P1, P2])
    return v


## K+ -----
t = np.arange(0, 34, 0.001)
holding_p = - 70
holding_t = 4.5
P1s = np.arange(-50, 50+1, 10)
P1_t = 29.5
P2 = -70

plt.figure(0)
for P1 in P1s:
    k = voltage_clamp(t, holding_p, holding_t, P1, P1_t, P2)
    plt.plot(t, k)
plt.title('K+')
plt.ylabel('Clamp Voltage (mV)')
plt.xlabel('Time (sec)')


## Na2+ SSA -----
holding_p = - 100
holding_t = 10
P1s = np.arange(-85, 20+1, 5)
P1_t = 120 + holding_t
P2 = -100
t = np.arange(0, P1_t+10, 0.001)

plt.figure(1)
for P1 in P1s:
    k = voltage_clamp(t, holding_p, holding_t, P1, P1_t, P2)
    plt.plot(t, k)
plt.title('Na2+ SSA')
plt.ylabel('Clamp Voltage (mV)')
plt.xlabel('Time (ms)')


## Na2+ SSI -----
holding_p = - 100
holding_t = 10
P1s = np.arange(-140, -65+1, 5)
P1_t = 500 + holding_t
P2 = -20
t = np.arange(0, P1_t+10, 0.001)

plt.figure(2)
for P1 in P1s:
    k = voltage_clamp(t, holding_p, holding_t, P1, P1_t, P2)
    plt.plot(t, k)
plt.title('Na2+ SSI')
plt.ylabel('Clamp Voltage (mV)')
plt.xlabel('Time (ms)')


## Na2+ Rec -----
holding_p = - 100
holding_t = 10
P1 = -20
P1_ts = np.arange((50+holding_t+1), (50+holding_t+88+1), 3)
P2 = -100
t = np.arange(0, max(P1_ts)+1, 0.001)

plt.figure(3)
for P1_t in P1_ts:
    k = voltage_clamp(t, holding_p, holding_t, P1, P1_t, P2)
    plt.plot(t, k)
plt.title('Na2+ Rec')
plt.ylabel('Clamp Voltage (mV)')
plt.xlabel('Time (ms)')
