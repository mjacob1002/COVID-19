# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 17:31:05 2020

@author: Mathew Jacob
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import random


class Locality:

    def __init__(self, S0, I0, R0, dt, days, beta, gamma):
        # amount of entries in the np array
        self.beta = beta
        self.gamma = gamma
        self.size = int(days / dt)
        self.S = np.zeros(self.size + 1)
        self.I = np.zeros(self.size + 1)
        self.R = np.zeros(self.size + 1)
        self.N = S0 + I0 + R0
        self.S[0], self.I[0], self.R[0] = S0, I0, R0

    """
    Returns:
        derivatives for the S, I, R respectively as a tuple
    """

    def deriv(self, spot: int) -> tuple:
        return -self.beta * self.S[spot] * self.I[spot], self.beta * self.S[spot] * self.I[spot] - self.gamma * self.I[
            spot], self.gamma * self.I[spot]

    """
    psim is a helper method to the commute method
    
    Returns:
        a different number depending on whether the simulator picks:
            Susceptible individual: returns 0
            Infected individual: returns 1
            Removed individual: returns 2
    """

    def psim(self, index) -> int:
        x = random.random()
        if 0 <= x < self.S[index] / self.N:
            return 0
        elif self.S[index] / self.N <= x < (self.S[index] + self.I[index]) / self.N:
            return 1
        else:
            return 2

    def update_SIR(self, index: int):
        prev = index - 1
        rates = deriv(prev)
        self.S[index] = self.S[prev] + dt * rates[0]
        self.I[index] = self.I[prev] + dt * rates[1]
        self.R[index] = self.R[prev] + dt * rates[2]

    """
    transition between localities self -> y per unit time
    adjust the instance variables accordingly
    this function will be run before every state evaluation
    
    """


def commute(x: Locality, y: Locality, prop: float, pos: int):
    # runs through the loop a certain amount of times
    for i in range(int(prop * x.N)):
        if x.psim(pos) == 0:
            y.S[pos] += 1
            x.S[pos] -= 1
        elif x.psim(pos) == 1:
            y.I[pos] += 1
            x.I[pos] -= 1
        else:
            y.R[pos] += 1
            x.R[pos] -= 1
    # adjust the N values after the commute
    x.N -= 1
    y.N += 1


def state_change(x: Locality, y: Locality, pos):
    # randomizes the proportion of people coming from one locality to another
    a = random.random() * .08 + random.random() * .04
    commute(x, y, a, pos)
    b = random.random() * .08 + random.random() * .04
    commute(x, y, b, pos)


# parameters being used

index = 1
dt = .1
days = 30
N = int(days / dt)
t_vec = np.linspace(0, days, num=N)
b_bergen = 1.5
g_bergen = 1
b_nyc = 1.5
g_nyc = 1

bergen = Locality(100, 1, 0, dt, days, b_bergen, g_bergen)
nyc = Locality(100, 1, 0, dt, days, b_nyc, g_nyc)

for i in range(1, N + 1):
    state_change(bergen, nyc)
    bergen.update_SIR(i)
    nyc.update_SIR(i)

SB, IB, RB = bergen.S, bergen.I, bergen.R
SN, IN, RN = nyc.S, nyc.I, nyc.R

plt.plot(t_vec, IB, 'ro')
plt.xlabel("Time")
plt.ylabel("Number of Infected")
