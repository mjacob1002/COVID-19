from SIRSD_age import SIRSD
from SIRSD_age import randEvent
import numpy as np
import matplotlib.pyplot as plt
import copy


class SIRSDV(SIRSD):
    # beta is effective transmissible rate
    # gamma is the probability of going from I -> R
    # kappa is probability of going from R -> S
    # eta is probability of going from S -> V
    # dt is the time differential
    # days is the number of days being simulated
    def __init__(self, N: int, S0: int, I0: int, R0: int, V0: int, beta: float, gamma: float, kappa: float,
                 eta: float, dt: float, days: int, timedelay = 0):
        assert S0 + I0 + R0 + V0 == N
        super(SIRSDV, self).__init__(N, S0, I0, R0, beta, gamma, kappa, dt, days)
        self.V0 = V0
        self.eta = eta
        self.V = np.zeros(self.numPoints)
        self.V[0] = V0
        self.timedelay = timedelay

    # initialize the sets to be used in the simulation
    # for now, assuming that R0, D0, and V0 are all zero. Will implement to include those other values
    # for future simulations
    def initialize(self):
        mu = 38.4
        sigma = 20
        Sset = set()
        Iset = set()
        Rset = set()
        Dset = set()
        Vset = set()
        # create the age of the first person and add it to the set
        for i in range(self.I0):
            p0 = np.random.normal(mu, sigma)
            Iset.add(p0)
        # create an array of ages for the remaining number of people
        agesSusceptible = np.random.normal(mu, sigma, self.S0)
        # loop through te array of ages and add them to the susceptible set
        # operating under the assumption that R0, D0 are both
        for person in agesSusceptible:
            Sset.add(person)
        assert len(Sset) == self.S0
        assert len(Iset) == self.I0
        return Sset, Iset, Rset, Dset, Vset

    def run(self):
        Sset, Iset, Rset, Dset, Vset = self.initialize()
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
                else:
                    # the if statement directly reflects a time delay of only when it is past day 30 to start vaccinating
                    # remove this if statement when time delay is not needed
                    if i * self.dt >= self.timedelay:
                        event = randEvent(self.eta)
                        # if the person does become vaccinated
                        if event:
                            Sset.remove(j)
                            Vset.add(j)
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
            self.S[i], self.I[i], self.R[i], self.D[i], self.V[i] = len(Sset), len(Iset), len(Rset), len(Dset), len(
                Vset)

    def plot(self):
        t = np.linspace(0, self.days, num=self.numPoints)
        plt.plot(t, self.S, label="Susceptibles", color='r')
        plt.plot(t, self.I, label="Infecteds", color='b')
        plt.plot(t, self.R, label="Recovereds", color='g')
        plt.plot(t, self.D, label="Deaths", color='m')
        plt.plot(t, self.V, label="Vaccinated", color='y')
        plt.xlabel("Number of Days")
        plt.ylabel("Number of People")
        plt.title("Epidemic Spread in a community, N=1000")
        plt.legend()
        plt.show()


N = 10000
S0 = 9950
I0 = 50
R0 = 0
V0 = 0
dt = .1
days = 120
beta = .8 * dt
gamma = .3 * dt
kappa = .053 * dt
eta = .05 * dt
test = SIRSDV(N=N, S0=S0, I0=I0, R0=R0, V0=V0, beta=beta, gamma=gamma, kappa=kappa, eta=eta, dt=dt, days=days, timedelay=90)
test.run()
test.plot()
