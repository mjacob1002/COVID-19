from SIRSD_periodic_spatial import SIRSDPeriodicSpatial, randEvent, randMovementAngle
import math
import numpy as np
from matplotlib import pyplot as plt
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


class Person:
    def __init__(self, x, y, age, ss, r, R, theta):
        # randomly generated x coordinate
        self.x = x
        # randomly generated y coordinate
        self.y = y
        # randomly generated age
        self.age = age
        # randomly generated super spreader boolean
        self.ss = ss
        # spreading radius
        self.r = r
        # randomly generated radius of movement
        self.R = R
        # randomly generated starting angle
        self.theta = theta
        self.currTheta = theta
        self.isIncluded = False
        # determined center from equation h = x - r*cos(theta) and k = y - r*sin(theta)
        self.h = self.x - R * math.cos(theta)
        self.k = self.y - R * math.sin(theta)


# eta is S -> V, zeta is leaving E, ioda is ICU->D, rho is ICU->R
class ICUVSpatial(SIRSDPeriodicSpatial):
    def __init__(self, N, S0, I0, R0, E0, V0, rstart, rstartsigma, days, gamma, kappa, eta, zeta, rho, phi, omega, ioda,
                 density,
                 planeSize, alpha=2):
        super(ICUVSpatial, self).__init__(N, S0, I0, R0, rstart, days, gamma, kappa, density, planeSize, 1,
                                          alpha=2)
        self.V, self.E = np.zeros(days + 1), np.zeros(days + 1)
        self.Lag = np.zeros(days + 1)
        self.ICU = np.zeros(days + 1)
        self.infectious = np.zeros(days + 1)
        self.combine = np.zeros(days + 1)
        self.additionalDeath = np.zeros(days + 1)
        self.V[0] = V0
        self.E[0] = E0
        self.Lag[0] = 0
        self.ICU[0] = 0
        self.infectious[0] = self.Lag[0] + self.I[0]
        self.Vcollect = []
        self.Lagcollect = []
        self.ICUcollect = []
        self.Ecollect = []
        self.rstartsigma = rstartsigma
        self.muAges = 45.8
        self.sigmaAges = 15
        self.muR = 12
        self.sigmaR = 5
        self.eta, self.zeta, self.rho, self.phi, self.omega, self.ioda = eta, zeta, rho, phi, omega, ioda

    def infect(self, p1, p2):
        r0 = p1.r
        r = ((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2) ** 0.5
        # r = math.sqrt(math.pow(p1.x-p2.x, 2) + math.pow(p1.y-p2.y, 2))
        if r >= r0:
            return False
        else:
            w = self.w0 * (1 - math.pow((r / r0), self.alpha))
            infectEvent = randEvent(w)
            return infectEvent

    # P(E->L|E->)
    def PICU(self, age: float):
        # bunch of sigmoid function stuff
        a = math.exp(-(age - 65))
        w0 = 1
        return w0 / (1 + a)

    def initialize(self):
        ages = np.random.normal(self.muAges, self.sigmaAges, self.N)
        coordinateX = np.random.random(self.N) * self.L
        coordinateY = np.random.random(self.N) * self.L
        movementRadii = np.random.normal(self.muR, self.sigmaR, self.N)
        startingAngles = np.random.random(self.N) * 2 * math.pi
        spreadingRadius = np.random.normal(self.rstart, self.rstartsigma, self.N)
        for i in range(self.N):
            # S copy
            p1 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), spreadingRadius[i],
                        movementRadii[i],
                        startingAngles[i])
            # I copy
            p2 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), spreadingRadius[i],
                        movementRadii[i],
                        startingAngles[i])
            # R copy
            p3 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), spreadingRadius[i],
                        movementRadii[i],
                        startingAngles[i])
            # D copy
            p4 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), spreadingRadius[i],
                        movementRadii[i],
                        startingAngles[i])
            # Ecopy
            p5 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), spreadingRadius[i],
                        movementRadii[i],
                        startingAngles[i])
            # Lcopy
            p6 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), spreadingRadius[i],
                        movementRadii[i],
                        startingAngles[i])
            # ICU copy
            p7 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), spreadingRadius[i],
                        movementRadii[i],
                        startingAngles[i])
            # Vaccine copy
            p8 = Person(coordinateX[i], coordinateY[i], ages[i], randEvent(self.density), spreadingRadius[i],
                        movementRadii[i],
                        startingAngles[i])
            # Susceptible collect
            self.Scollect.append(p1)
            # Exposed collect
            self.Ecollect.append(p5)
            # Infected collected
            self.Icollect.append(p2)
            # Lag compartment collected
            self.Lagcollect.append(p6)
            # Removed compartment collected
            self.Rcollect.append(p3)
            # ICU compartment collect
            self.ICUcollect.append(p7)
            # Death compartment collect
            self.Dcollect.append(p4)
            # Vaccinated compartment collect
            self.Vcollect.append(p8)
        for i in range(self.N):
            # print("Iterator:", i, " I[0]: ", self.I[0], "I[i].isIncluded: ", self.Icollect[i].isIncluded, "Condition:", i<self.I[0])
            if i < self.I[0]:
                self.Icollect[i].isIncluded = True
            elif i < self.I[0] + self.S[0]:
                self.Scollect[i].isIncluded = True
            elif i < self.I[0] + self.S[0] + self.V[0]:
                self.Vcollect[i].isIncluded = True
            else:
                self.Rcollect[i].isIncluded = True
            # print("Icollect[i].isIncluded:", self.Icollect[i].isIncluded)

    # runs in N^2 time
    # does the state changes from S to I
    # returns a set of the indicies which'll into I
    def StoE(self):
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
            # same process for the lag group, as they're also infectious
            for j in self.Lagcollect:
                if not j.isIncluded:
                    continue
                event = self.infect(j, i)
                if not event:
                    continue
                i.isIncluded = False
                transfers.add(count)

        # return the set so that the transfers can be done at the end
        # this is because at one state change, you don't people who just turned into I from S being capable of
        # immediately going from E to I/L
        return transfers

    # run state changes for those leaving the E compartment
    def leaveE(self):
        # keeps track of those going from I - > R
        transfersL = set()
        transfersI = set()
        for count, person in enumerate(self.Ecollect):
            # if the person isn't an infected person
            if not person.isIncluded:
                continue
            # if the person is infected
            # test for leaving E
            w = randEvent(self.zeta)
            leaveE = randEvent(w)
            # if the person leaves E, test for whether they go to L or I
            if not leaveE:
                continue
            # remove them from infected
            person.isIncluded = False
            # test whether they'll go to L or I using PICU
            # get the probability
            w = self.PICU(person.age)
            # generate an event
            eventICU = randEvent(w)
            # if the person is determined to go to the ICU
            if eventICU:
                transfersL.add(count)
            else:
                transfersI.add(count)
        # those who wil go from I -> R, but not done immediately because state change will be done on current R
        # and those who just transfered shouldn't be eligible to go immediately back to S
        return transfersL, transfersI

    # run state changes for those leaving L -> ICU
    def LtoICU(self):
        transfersICU = set()
        for count, person in enumerate(self.Lagcollect):
            # if the person currently isn't in lag compartment, move onto the next
            if not person.isIncluded:
                continue
            # person is in lag compartment, so generate an event
            event = randEvent(self.phi)
            # if the person doesn't transfer to ICU, then don't do anything and check next person
            if not event:
                continue
            # the person does transfer to ICU
            # remove them from the Lcollect
            self.Lagcollect[count].isIncluded = False
            # add the index to the set to be changed at the end
            transfersICU.add(count)
        return transfersICU

    # run state changes for those leaving I -> R/D
    # pos is the current day the state change is for
    def leaveI(self, pos):
        # set containing going to D
        transfersD = set()
        # set contianing those leaving R
        transfersR = set()
        for count, person in enumerate(self.Icollect):
            # check to make sure the person is included in the Icollect
            if not person.isIncluded:
                continue
            # check the death event first
            death = randEvent(self.omega)
            # if the person dies
            if death:
                # change the person from being included in I to D
                self.Icollect[count].isIncluded = False
                transfersD.add(count)
                # move on to the next infected person
                continue
            # if the person hasn't died, check to see whether they've recovered by creating an evnet
            recovery = randEvent(self.gamma)
            # if the person hasn't recovered, nothing to do. Just continue on to the next infected
            if not recovery:
                continue
            # if the person has recovered, remove the person from Icollect and add index to the set transferR
            self.Icollect[count].isIncluded = False
            transfersR.add(count)
        return transfersD, transfersR

    # pos is the current index value of the state change i.e day 1, day 5, etc
    def leaveICU(self, pos):
        transfersD = set()
        transfersR = set()
        for count, person in enumerate(self.ICUcollect):
            # determine if the person is in the ICU collection
            if not person.isIncluded:
                # move on to the next person, as current person isn't in ICU
                continue
            # determine whether the person dies
            death = randEvent(self.ioda)
            # if the person dies
            if death:
                # set isIncluded flag to False in ICUcollect and True in Death Collection
                self.ICUcollect[count].isIncluded = False
                transfersD.add(count)
                # move on to the next person in the ICU
                continue
            # if the person hasn't died, determine if they've recovered
            recovery = randEvent(self.rho)
            # if ther person hasn't recovered, run checks on the next ICU patient
            if not recovery:
                continue
            # if the person has recovered,
            self.ICUcollect[count].isIncluded = False
            transfersR.add(count)
        return transfersD, transfersR

    # run S -> V
    # pos is the current day the simulation is taking place on
    def StoV(self, pos):
        # keeps track of the transfers
        transfers = set()
        for count, person in enumerate(self.Scollect):
            # check to see if person is included in the Scollection
            # print("Evaluating person.isIncluded:", person.isIncluded)
            if not person.isIncluded:
                # print("False")
                continue
            # if the person is in susceptible, generate a vaccination event
            vaccination = randEvent(self.eta)
            # if the person isn't vaccinated, continue on to the next person
            if not vaccination:
                continue
            # if the person becomes vaccinated
            # change the collections, where false in susceptible and true in vaccinated
            self.Scollect[count].isIncluded = False
            self.Vcollect[count].isIncluded = True
            # change the amount kept track of in the array
            transfers.add(count)
        return transfers

    # run a verification to make sure that there is a single isIncluded for every person
    def assertionCheck(self):

        for i in range(self.N):
            # the current number of trues in a given row
            trueCount = 0
            # run checks to get a count on the number of trues in a single row
            if self.Scollect[i].isIncluded:
                trueCount += 1
            if self.Ecollect[i].isIncluded:
                trueCount += 1
            if self.Icollect[i].isIncluded:
                trueCount += 1
            if self.Lagcollect[i].isIncluded:
                trueCount += 1
            if self.ICUcollect[i].isIncluded:
                trueCount += 1
            if self.Rcollect[i].isIncluded:
                trueCount += 1
            if self.Dcollect[i].isIncluded:
                trueCount += 1
            if self.Vcollect[i].isIncluded:
                trueCount += 1
            assert trueCount == 1

    def move(self):
        # for all of the people in the population
        for i in range(self.N):
            # generate a random movement angle for the person at
            angle = randMovementAngle()
            # if the person isn't dead, then adjust his coordinates.
            # if the person is dead, he won't move nor leave the D state, so nothing needs to be done
            if not self.Dcollect[i].isIncluded:
                self.adjustTheta(self.Scollect[i], angle)
                self.adjustTheta(self.Ecollect[i], angle)
                self.adjustTheta(self.Icollect[i], angle)
                self.adjustTheta(self.Lagcollect[i], angle)
                self.adjustTheta(self.Rcollect[i], angle)
                self.adjustTheta(self.Vcollect[i], angle)
                self.adjustTheta(self.Dcollect[i], angle)

    def stateChanger(self, transfersStoE, transfersEtoL, transfersEtoI, transfersItoD, transfersItoR, transfersLtoICU,
                     transfersICUtoD, transfersICUtoR, transfersRtoS, transfersStoV):
        # run the isIncluded changes
        for index in transfersStoE:
            self.Scollect[index].isIncluded = False
            self.Ecollect[index].isIncluded = True
        for index in transfersEtoL:
            self.Ecollect[index].isIncluded = False
            self.Lagcollect[index].isIncluded = True
        for index in transfersEtoI:
            self.Ecollect[index].isIncluded = False
            self.Icollect[index].isIncluded = True
        for index in transfersItoR:
            self.Icollect[index].isIncluded = False
            self.Rcollect[index].isIncluded = True
        for index in transfersLtoICU:
            self.Lagcollect[index].isIncluded = False
            self.ICUcollect[index].isIncluded = True
        for index in transfersICUtoR:
            self.ICUcollect[index].isIncluded = False
            self.Rcollect[index].isIncluded = True
        for index in transfersRtoS:
            self.Rcollect[index].isIncluded = False
            self.Scollect[index].isIncluded = True
        for index in transfersItoD:
            self.Dcollect[index].isIncluded = True
        for index in transfersICUtoD:
            self.Dcollect[index].isIncluded = True
        for index in transfersStoV:
            self.Vcollect[index].isIncluded = True

    @profile
    def run(self):
        self.initialize()
        for i in range(1, days + 1):
            # run the state changes for S to E
            transfersStoE = self.StoE()
            transfersEtoL, transfersEtoI = self.leaveE()
            # run state change from I to R/D
            transfersItoD, transfersItoR = self.leaveI(i)
            # run state change from L to ICU
            transfersLtoICU = self.LtoICU()
            # run state change from ICU to R/D
            transfersICUtoD, transfersICUtoR = self.leaveICU(i)
            # run state change for R to S
            transfersRtoS = self.RtoS()
            # run the state changes for S to V
            transfersStoV = self.StoV(i)

            self.stateChanger(transfersStoE, transfersEtoL, transfersEtoI, transfersItoD, transfersItoR,
                              transfersLtoICU, transfersICUtoD, transfersICUtoR, transfersRtoS, transfersStoV)

            # make sure that everything functions properly
            self.assertionCheck()
            # adjust values; don't need to factor in absorption states as they're done during state change
            # S adjustment doesn't include V because V is absorption state
            self.S[i] += (self.S[i - 1] + len(transfersRtoS) - len(transfersStoE) - len(transfersStoV))
            self.E[i] += (self.E[i - 1] + len(transfersStoE) - (len(transfersEtoI) + len(transfersEtoL)))
            self.I[i] += (self.I[i - 1] + len(transfersEtoI) - len(transfersItoR) - len(transfersItoD))
            self.Lag[i] += (self.Lag[i - 1] + len(transfersEtoL) - len(transfersLtoICU))
            self.ICU[i] += (self.ICU[i - 1] + len(transfersLtoICU) - len(transfersICUtoR) - len(transfersICUtoD))
            self.R[i] += (self.R[i - 1] + len(transfersICUtoR) + len(transfersItoR) - len(transfersRtoS))
            self.V[i] += self.V[i - 1] + len(transfersStoV)
            self.D[i] += self.D[i - 1] + len(transfersItoD) + len(transfersICUtoD)
            self.infectious[i] = self.I[i] + self.Lag[i]
            self.combine[i] = self.combine[i] + self.ICU[i]
            # move everyone after state change is complete
            self.move()
            # make sure that the total number of people in system remains constant
            assert self.S[i] + self.E[i] + self.I[i] + self.Lag[i] + self.R[i] + self.ICU[i] + self.D[i] + self.V[
                i] == self.N

    def plot(self):
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex='all')
        t = np.linspace(0, self.days, num=self.days + 1)
        # plt.plot(t, self.S, label="Susceptibles", color='r')
        ax1.plot(t, self.infectious, label="Infectious", color='b')
        ax1.set_ylabel("Active Cases")
        ax1.set_title("Epidemic Model of N = 2000")
        ax1.legend()
        # plt.plot(t, self.R, label="Recovereds", color='g')
        ax2.plot(t, self.D, label="Deaths", color='m')
        ax2.set_ylabel("Cumulative Deaths")
        ax2.legend()
        # plt.plot(t, self.E, label="Exposed", color='r')
        ax3.plot(t, self.ICU, label="ICU")
        ax3.set_xlabel("Days")
        ax3.set_ylabel("Active Hospitalizations")
        ax3.legend()
        # plt.plot(t, self.V, label="Vaccinated")
        plt.show()

    def testing(self):
        for i in range(1, days + 1):
            self.StoV(i)

    def maximizer(self):
        maximum = [self.infectious[0], self.ICU[0]]
        index = np.zeros(2)
        for i in range(self.days + 1):
            curr1 = self.infectious[i]
            maximum[0] = max(maximum[0], curr1)
            if maximum[0] == curr1:
                index[0] = i
            curr2 = self.ICU[i]
            maximum[1] = max(maximum[1], curr2)
            if maximum[1] == curr2:
                index[1] = i

        return index


N = 1000
S0 = 999
I0 = 1
R0 = 0
D0 = 0
V0 = N - S0 - I0
rstart = 2
sigmaR = .5
L = 50  # 125
gamma = .3
kappa = 0.1
omega = 0.001
zeta = .2
phi = .3
rho = .3
ioda = .2
eta = 0.007
density = .33
days = 120

test = ICUVSpatial(N, S0, I0, R0, 0, V0, rstart, sigmaR, days, gamma, kappa, eta, zeta, rho, phi, omega, ioda, density,
                   L)
test.run()
print(test.maximizer())
test.plot()
#print(test.D)


# test.run()
# test.plot()
