# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 17:31:05 2020

@author: Mathew Jacob
"""

import numpy as np
from matplotlib import pyplot as plt
import random


class Locality:
    """Creates an object of a locality, such as New York City or Bergen County.
    Attributes:
        beta: the effective transmission of the locality
        gamma: the average recovery rate of those with COVID-19 in the locality
        S0: the starting amount of susceptible individuals at start of simulation
        I0: the starting amount of infected individuals at start of simulation
        R0: the starting amount of removed/recovered individuals at start of simulation
        dt: the differential used for Euler's method.
        days: amount of days that are being simulated
        size: the total amount of points being calculated
        N: the total amount of people in the locality at any given time
    """
    def __init__(self, beta, gamma, S0, I0, R0, dt, days):
        """
        The constructor for the Locality class.
        :param beta: the effective transmission of the locality
        :param gamma: the average recovery rate of those with COVID-19 in the locality
        :param S0: he starting amount of susceptible individuals at start of simulation
        :param I0: the starting amount of infected individuals at start of simulation
        :param R0: the starting amount of removed/recovered individuals at start of simulation
        :param dt: the differential used for Euler's method.
        :param days: amount of days that are being simulated
        """
        self.size = int(days / dt)
        self.dt = dt
        self.days = days
        self.S = np.zeros(self.size + 1)
        self.I = np.zeros(self.size + 1)
        self.R = np.zeros(self.size + 1)
        self.beta = beta
        self.gamma = gamma
        self.N = S0 + I0 + R0
        self.S[0], self.I[0], self.R[0] = S0, I0, R0

    def f1(self, pos: int):
        """Computes the derivative at a given time
        Paramters:
        pos(int): the element in S, I, R that the derivative is calculated for

        Returns:
            the derivatives of S(t), I(t), and R(t) as a tuple
        """
        x = -1 * self.beta * self.S[pos] * self.I[pos] / self.N
        y = self.beta * self.S[pos] * self.I[pos] / self.N - self.gamma * self.I[pos]
        z = self.gamma * self.I[pos]
        return x, y, z

    def psim(self, index: int) -> int:
        """Helper method of commute()that simulates the state of people that commute from one locality to another
        Parameter:
            index(int) : the current 'time' that it is aka the latest index in the S, I, R arrays
        Returns:
            0 if the simulation determines that the person is susceptible, 1 for infected, and 2 for a recovered person
        """
        # p is probability that the person commuting is from Susceptible
        p = self.S[index] / self.N
        # q is probability that the person commuting is from Infected'
        q = self.I[index] / self.N
        a = random.random()
        if 0 <= a < p:
            return 0
        elif p <= a < p + q:
            return 1
        else:
            return 2

    def update_sir(self, pos: int):
        """
        Uses Euler's Method to create an accurate approximation of the spread of COVID
        :param pos: the current cell in S, I, and R that is being simulate/ calculated
        :return:
        """
        change: tuple
        change = self.f1(pos - 1)
        self.S[pos] = self.S[pos - 1] + self.dt * change[0]
        self.I[pos] = self.I[pos - 1] + self.dt * change[1]
        self.R[pos] = self.R[pos - 1] + self.dt * change[2]


def commute(x: Locality, y: Locality, index, prop: float):
    """
    Simulates the commute of people from locality x to locality y
    :param x: The first locality. People leave from locality x to locality y
    :param y: The second locality. People come to locality y from locality x
    :param prop: the proportion of individuals that go from x -> y
    :param index: The current index of the S, I, R values that are being changed.
    :return: nothing; this is a void function
    """

    # simulates the commuting process for the given amount of people that are supposed to commute
    for i in range(int(prop * x.N)):
        if x.psim(index) == 0:
            y.S[index] += 1
            x.S[index] -= 1
        elif x.psim(index) == 1:
            y.I[index] += 1
            x.I[index] -= 1
        else:
            y.R[index] += 1
            x.R[index] -= 1
        # change the total population of a locality accordingly
        y.N += 1
        x.N -= 1


def state_change(x: Locality, y: Locality, index: int, prop1: float, prop2: float):
    """
    Implements the commute and update_sir methods to actual simulate commute and spread of disease in locality.
    :param x: one of two localities used in this method
    :param y: second of two localities used in this method
    :param prop1: the proportion of people who go from x -> y
    :param prop2: the proportion of people who go from y -> x
    :param index: the current cell of S, I, R arrays that are being calculated using update_sir(index)
    :return:
    """
    # copying the existing value to be reapplied after the commutes and evaluation
    pos = index - 1
    sx, ix, rx = x.S[pos], x.I[pos], x.R[pos]
    sy, iy, ry = y.S[pos], y.I[pos], y.R[pos]
    # does the commute between localities and changes the values that were copied over
    commute(x, y, pos, prop1)
    commute(y, x, pos, prop2)
    # updates the next cell using the new S, I, R values
    x.update_sir(index)
    y.update_sir(index)
    # reapplies the copied values to where they should be
    x.S[pos], x.I[pos], x.R[pos] = sx, ix, rx
    y.S[pos], y.I[pos], y.R[pos] = sy, iy, ry


beta1 = 1.5
gamma1 = 3/15
S01, I01, R01 = 999, 1, 0

beta2 = 2
gamma2 = 3/15
S02, I02, R02 = 50, 1, 0

day = 31
dt = .1
num_times = int(day / dt)
t = np.linspace(0, day, num=num_times+1)

nyc = Locality(beta=beta1, gamma=gamma1, S0=S01, I0=I01, R0=R01, dt=dt, days=day)
bergen = Locality(beta=beta2, gamma=gamma2, S0=S02, I0=I02, R0=R02, dt=dt, days=day)
for i in range(1, num_times + 1):
    state_change(bergen, nyc, i,prop1=.5, prop2 = .1 )
    print("I: " , nyc.I[i])
plt.plot(t, nyc.I, 'r-', label="NYC Infected")
plt.plot(t, bergen.I, 'b-', label="Bergen County Infected")
plt.xlabel("Days")
plt.ylabel("Number of Infected")
plt.title("Locality Model: Simulation 4")
plt.legend()
plt.show()
