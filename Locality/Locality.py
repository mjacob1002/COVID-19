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
        #amount of entries in the np array
        self.beta = beta
        self.gamma = gamma
        self.size = int(days/dt)
        self.S = np.zeros(self.size+1)
        self.I = np.zeros(self.size+1)
        self.R = np.zeros(self.size+1)
        self.N = S0 + I0 + R0
        self.t[0], self.S[0], self.I[0], self.R[0] = 0, S0, I0, R0
        
    """
    Returns:
        derivatives for the S, I, R respectively as a tuple
    """
    def deriv(self, spot : int) -> tuple:
        return (-self.beta*self.S[spot] * self.I[spot], self.beta * self.S[spot] * self.I[spot] - self.gamma* self.I[spot], self.gamma* self.I[spot])
   
    """
    psim is a helper method to the commute method
    
    Returns:
        a different number depending on whether the simulator picks:
            Susceptible individual: returns 0
            Infected individual: returns 1
            Removed individual: returns 2
    """
    def psim(self) -> int:
        x = random.random()
        if 0 <= x < (self.S)/self.N:
            return 0
        elif self.S/self.N <= x < (self.S+self.I)/self.N:
            return 1
        else:
            return 2

    def update_SIR(self, index: int):
        prev = index -1
        rates = deriv(prev)
        self.S[index] = self.S[prev] + dt * rates[0]
        self.I[index] = self.I[prev] + dt * rates[1]
        self.R[index] = self.R[prev] + dt * rates[2]
    
    

    """
    transition between localities self -> y per unit time
    adjust the instance variables accordingly
    this function will be run before every state evaluation
    
    """
def commute(x: Locality, y: Locality, prop: float):
        # runs through the loop a certain amount of times
        for i in range(int(prop*x.N)):
            if(x.psim() == 0):
                y.S +=1
                x.S -=1
            elif x.psim() == 1:
                y.I +=1
                x.I -=1
            else:
                y.R +=1
                x.R -=1
        #adjust the N values after the commute
        x.N-=1
        y.N +=1 
# uses the commute method to commute people from x -> y as well as y -> x 
def state_change(x: Locality, y: Locality):
    # randomizes the proportion of people coming from one locality to another
    a = random.random() * .08 + random.random()*.04
    commute(x, y, a)
    b = random.random() * .08 + random.random()*.04
    commute(x, y, b)

       
#parameters being used 
index = 1
dt = .1
days = 30
t_vec = np.zeros(days+1)
b_bergen = 1.5
g_bergen = 1
b_nyc = 1.5
g_nyc = 1

Bergen = Locality(100,1,0, dt, days, b_bergen, g_bergen)
NYC = Locality(100,1,0,dt, days, b_nyc, g_nyc)


for i in range(1, days+1):
    t_vec[i] = i
