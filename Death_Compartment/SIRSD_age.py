import numpy as np
import math
import copy
from matplotlib import pyplot as plt


def randEvent(p: float):
    if p <= 0:
        return False
    x = np.random.random()
    if 0 < x <= p:
        return True
    else:
        return False


class SIRSD:
    # beta is effective transmissible rate
    # gamma is the probability of going from I -> R
    # kappa is probability of going from R -> S
    # dt is the time differential
    # days is the number of days being simulated
    def __init__(self, N: int, S0: int, I0: int, R0: int, beta: float, gamma: float, kappa: float
                 , dt: float, days: int):
        self.N = N
        self.S0 = S0
        self.I0 = I0
        self.R0 = R0
        self.D0 = 0
        self.beta = beta
        self.gamma = gamma
        self.kappa = kappa
        self.dt = dt
        self.days = days
        self.numPoints = int(days / dt + 1)
        self.S, self.I, self.R, self.D = np.zeros(self.numPoints), np.zeros(self.numPoints), np.zeros(self.numPoints), np.zeros(
            self.numPoints)
        self.S[0], self.I[0], self.R[0], self.D[0] = S0, I0, R0, self.D0

    def PDeath(self, age: float):
        # bunch of sigmoid function stuff
        a = math.exp(-(age - 65))
        w0 = .73
        return w0 / (1 + a)

    def infect(self, i, pop):
        return self.beta * i / pop

    # initialize the sets to be used in the simulation
    def initialize(self):
        mu = 38.4
        sigma = 20
        Sset = set()
        Iset = set()
        Rset = set()
        Dset = set()
        # create the age of the first person and add it to the set
        for i in range(self.I0):
            p0 = np.random.normal(mu, sigma)
            Iset.add(p0)
        # create an array of ages for the remaining number of people
        ages = np.random.normal(mu, sigma, self.N - self.I0)
        # loop through te array of ages and add them to the susceptible set
        # operating under the assumption that R0, D0 are both
        for person in ages:
            Sset.add(person)
        assert len(Sset) == self.S0
        assert len(Iset) == self.I0
        return Sset, Iset, Rset, Dset

    def run(self):
        Sset, Iset, Rset, Dset = self.initialize()
        for i in range(1, self.numPoints):
            # make copies of the sets that will have removals
            Scopy = copy.deepcopy(Sset)
            Icopy = copy.deepcopy(Iset)
            Rcopy = copy.deepcopy(Rset)
            mu = self.infect(self.I[i - 1], self.N)
            for j in Scopy:
                # generate a random event determining whether someone becomes infected
                event = randEvent(mu)
                # if the person does become infected
                if event:
                    Sset.remove(j)
                    Iset.add(j)
            for j in Icopy:
                # determine if someone dies
                dprob = self.PDeath(j)
                event = randEvent(dprob)
                # if someone dies, modify both the sets and the total population size
                if event:
                    Iset.remove(j)
                    Dset.add(j)
                    self.N -= 1
                # determine if someone recovers if they don't die
                else:
                    event2 = randEvent(self.gamma)
                    if event2:
                        Iset.remove(j)
                        Rset.add(j)
            for j in Rcopy:
                event = randEvent(self.kappa)
                if event:
                    Rset.remove(j)
                    Sset.add(j)
            self.S[i], self.I[i], self.R[i], self.D[i] = len(Sset), len(Iset), len(Rset), len(Dset)

    def plot(self):
        t = np.linspace(0, self.days, num=self.numPoints)
        plt.plot(t, self.S, label="Susceptibles", color='r')
        plt.plot(t, self.I, label="Infecteds", color='b')
        plt.plot(t, self.R, label="Recovereds", color='g')
        plt.plot(t, self.D, label="Deaths", color='m')
        plt.xlabel("Number of Days")
        plt.ylabel("Number of People")
        plt.title("Epidemic Spread in a community, N=1000")
        plt.legend()
        plt.legend()
        plt.plot()
        plt.show()


N = 1000
S0 = 950
I0 = 50
R0 = 0
dt = .1
days = 120
beta = .8 * dt
gamma = .3 * dt
kappa = .0053 * dt

test = SIRSD(N=N, S0=S0, I0=I0, R0=R0, beta=beta, gamma=gamma, kappa=kappa, dt=dt, days=days)
#test.run()
#test.plot()
