import math

import numpy as np
from multipledispatch import dispatch
import copy
from matplotlib import pyplot as plt
import time
import cProfile, pstats, io


def profile(fnc):
    """A decorator that uses cProfile to profile a function"""

    def inner(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner


def randEvent(p: float):
    if p <= 0:
        return False
    x = np.random.random()
    if 0 < x <= p:
        return True
    else:
        return False


def randMovementAngle():
    delta_theta = np.random.uniform(-2 * math.pi, 2 * math.pi)
    return delta_theta


class Person:
    def __init__(self, x, y, age, ss, R, theta):
        # randomly generated x coordinate
        self.x = x
        # randomly generated y coordinate
        self.y = y
        # randomly generated age
        self.age = age
        # randomly generated super spreader boolean
        self.ss = ss
        # randomly generated radius of movement
        self.R = R
        # randomly generated starting angle
        self.theta = theta
        self.currTheta = theta
        self.isIncluded = False
        # determined center from equation h = x - r*cos(theta) and k = y - r*sin(theta)
        self.h = self.x - R * math.cos(theta)
        self.k = self.y - R * math.sin(theta)

    def __repr__(self):
        return "Person(%s,%s)" % (self.x, self.y)

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        if isinstance(other, Person):
            return (self.x == other.x) and (self.y == other.y)
        else:
            return False


# trying to make this spatial model faster by not using deepcopy a lot, as that takes most of the time in the other version
# runs in tN^2 time, where t is the number of days simulated, and N is the population size

class SIRSDPeriodicSpatial:
    def __init__(self, N, S0, I0, R0, rstart, days, gamma, kappa, density, L, hubConstant=6 ** .5, alpha=2):
        # size of the population
        self.N = N
        # radius of infection
        self.rstart = rstart
        # number of days
        self.days = days
        # probability of going from I -> R
        self.gamma = gamma
        # probability of going from R -> S
        self.kappa = kappa
        # probability of becoming a super spreader
        self.density = density
        # the size of the plane
        self.L = L
        # arrays containing the number of susceptibles, infecteds, removeds ,and dead people
        self.S, self.I, self.R, self.D = np.zeros(days + 1), np.zeros(days + 1), np.zeros(days + 1), np.zeros(days + 1)
        self.S[0], self.I[0], self.R[0], self.D[0] = S0, I0, R0, 0
        # lists that will contain the Person objects; will act as a replacement to the sets that were previously used
        self.Scollect, self.Icollect, self.Rcollect, self.Dcollect = [], [], [], []
        self.w0 = 1
        self.hubConstant = hubConstant
        self.alpha = alpha
        # temporary variable to confirm that people are moving
        self.temp = []

    # @dispatch(Person, Person)
    def infect(self, p1, p2):
        r0: float
        if p1.ss:
            r0 = self.rstart * self.hubConstant
        else:
            r0 = self.rstart
        r = ((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2) ** 0.5
        # r = math.sqrt(math.pow(p1.x-p2.x, 2) + math.pow(p1.y-p2.y, 2))
        if r >= r0:
            return False
        else:
            w = math.pow(self.w0 * (1 - r / r0), self.alpha)
            infectEvent = randEvent(w)
            return infectEvent

    def initialize(self):
        mu = 45.4
        sigma = 15
        muR = 12
        sigmaR = 10
        ages = np.random.normal(mu, sigma, self.N)
        coordinateX = np.random.random(self.N) * self.L
        coordinateY = np.random.random(self.N) * self.L
        movementRadii = np.random.normal(muR, sigmaR, self.N)
        startingAngles = np.random.random(self.N) * 2 * math.pi
        for i in range(self.N):
            p1 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), movementRadii[i],
                        startingAngles[i])
            p2 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), movementRadii[i],
                        startingAngles[i])
            p3 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), movementRadii[i],
                        startingAngles[i])
            p4 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), movementRadii[i],
                        startingAngles[i])
            self.Scollect.append(p1)
            self.Icollect.append(p2)
            self.Rcollect.append(p3)
            self.Dcollect.append(p4)
        for i in range(self.N):
            # print("Iterator:", i, " I[0]: ", self.I[0], "I[i].isIncluded: ", self.Icollect[i].isIncluded, "Condition:", i<self.I[0])
            if i < self.I[0]:
                self.Icollect[i].isIncluded = True
            elif i < self.I[0] + self.S[0]:
                self.Scollect[i].isIncluded = True
            else:
                self.Rcollect[i].isIncluded = True
            # print("Icollect[i].isIncluded:", self.Icollect[i].isIncluded)

    def PDeath(self, age: float):
        # bunch of sigmoid function stuff
        a = math.exp(-(age - 65))
        w0 = .73
        return w0 / (1 + a)

    def adjustTheta(self, person: Person, delta: float):
        # adjust the current theta of the person object
        person.currTheta += delta
        # adjust the theoretical x values
        x = person.h + person.R * math.cos(person.currTheta)
        # perform boundary checks to ensure that x is between 0 <= x <= L
        if x > self.L:
            x = self.L
        elif x < 0:
            x = 0
        # adjust the y coordinate of the person
        y = person.k + person.R * math.sin(person.currTheta)
        # perform boundary checks on y to ensure that y is in [0,L]
        if y > self.L:
            y = self.L
        elif y < 0:
            y = 0
        person.x = x
        person.y = y

    def move(self):
        for i in range(self.N):
            # create the random angle by which the person has moved
            angle = randMovementAngle()
            # adjust the x and y coordinates according to the angle they've moved
            self.adjustTheta(self.Scollect[i], angle)
            self.adjustTheta(self.Icollect[i], angle)
            self.adjustTheta(self.Rcollect[i], angle)
            # only adjust the person's copy in the Dcollection array if they're not included to ensure dead people
            # aren't moving
            if not self.Dcollect[i].isIncluded:
                self.adjustTheta(self.Dcollect[i], angle)

    # runs in N^2 time
    # does the state changes from S to I
    # returns a set of the indicies which'll into I
    def StoI(self):
        # set to keep track of the transfers; keeps track of the indices that are to be transferred
        transfers = set()
        for count, i in enumerate(self.Scollect):
            # if the particular object isn't included in susceptible bin, skip over it
            if not i.isIncluded:
                continue
            for j in self.Icollect:
                # print("StoI:", type(j))
                # if the object isn't an infected person, don't bother doing any calculation
                if not j.isIncluded:
                    continue
                else:
                    # calculate the infection probability and determine the event
                    event = self.infect(j, i)
                    # if person is infected, turn the isIncluded off in S and keep track of index "count" in the set
                    if event:
                        # remove it from the Susceptible array
                        i.isIncluded = False
                        # put it in the transfers array to be converted to infected in the end
                        transfers.add(count)
                        # break the loop because already infected, no need to run calculations with other infected
                        # people
                        break
        # return the set so that the transfers can be done at the end
        # this is because at one state change, you don't people who just turned into I from S being capable of
        # immediately going from I to R/D
        return transfers

    # pass in the current day for the state change
    def leaveI(self, pos):
        # keeps track of those going from I - > R
        transfers = set()
        for count, person in enumerate(self.Icollect):
            # if the person isn't an infected person
            if not person.isIncluded:
                continue
            # if the person is infected
            else:
                # test for death event
                w = self.PDeath(person.age)
                death = randEvent(w)
                # person has died; immediately switch them from I -> R because D is absorption state
                if death:
                    # remove them from infected
                    person.isIncluded = False
                    # add them to the D array
                    self.Dcollect[count].isIncluded = True
                    # adjust the number of infected, the number of dead people
                    self.I[pos] -= 1
                    self.D[pos] += 1
                else:
                    recovery = randEvent(self.gamma)
                    if recovery:
                        person.isIncluded = False
                        transfers.add(count)

        # those who wil go from I -> R, but not done immediately because state change will be done on current R
        # and those who just transfered shouldn't be eligible to go immediately back to S
        return transfers

    def RtoS(self):
        transfers = set()
        for count, person in enumerate(self.Rcollect):
            if not person.isIncluded:
                continue
            else:
                resusceptible = randEvent(self.kappa)
                if resusceptible:
                    transfers.add(count)
                    self.Rcollect[count].isIncluded = False
        return transfers

    @profile
    def run(self):
        self.initialize()
        for i in range(1, self.days + 1):
            self.move()
            # temp function that should verify that people move; delete once confirmed
            self.temp.append((self.Scollect[0].x, self.Scollect[0].y))
            StoI = self.StoI()
            ItoR = self.leaveI(i)
            RtoS = self.RtoS()
            for index in StoI:
                self.Icollect[index].isIncluded = True
            for index in ItoR:
                self.Rcollect[index].isIncluded = True
            for index in RtoS:
                self.Scollect[index].isIncluded = True
            self.S[i] = self.S[i - 1] + self.S[i] - len(StoI) + len(RtoS)
            self.I[i] = self.I[i - 1] + self.I[i] - len(ItoR) + len(StoI)
            self.R[i] = self.R[i - 1] + self.R[i] - len(RtoS) + len(ItoR)
            self.D[i] += self.D[i - 1]

    def printType(self):
        for i in self.Icollect:
            print(type(i))

    def printIsIncluded(self):
        for i in self.Icollect:
            print(i.isIncluded)

    def plot(self):
        t = np.linspace(0, self.days, num=self.days + 1)
        plt.plot(t, self.S, label="Susceptibles", color='r')
        plt.plot(t, self.I, label="Infecteds", color='b')
        plt.plot(t, self.R, label="Recovereds", color='g')
        plt.plot(t, self.D, label="Deaths", color='m')
        plt.xlabel("Number of Days")
        plt.ylabel("Number of People")
        plt.title("Epidemic Spread in a community, N=1000")
        plt.legend()
        plt.show()

    def tempFunction(self):
        print(self.temp)


#N = 10000
#S0 = 9950
#I0 = N - S0
#R0 = 0
#D0 = 0
#rstart = .5
#L = 100
#gamma = .45
#kappa = 0.3
#density = .33
#days = 31

#test = SIRSDPeriodicSpatial(N, S0, I0, R0, rstart, days, gamma, kappa, density, L)
#test.run()
# test.printIsIncluded()
# test.run()
#test.plot()
# print(test.R)
# print(test.D)
#test.tempFunction()
