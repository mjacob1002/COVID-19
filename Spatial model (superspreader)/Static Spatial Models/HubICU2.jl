using Random
using Plots

# Parameters
# time simulated
t = 80
# the exponent for the infection probability function
α = 2
# the normal radius of a non-super spreader
rstart = 6.0
# total population
N = 100000
# size of one side of the plane
sizeL = 250
# starting probability when two people are on top of each other
w0 = 1.0
# density of super spreaders
λ = .3
# probability of exposed to infectious(zeta)
ζ = .3
# probability of, once infected, going to L compartment
ρ = .03
# probability of going from I to R
γ = .3
# probability of going from L to ICU
η = .3
# probability of going from ICU to R
κ = .3
# boolean to determine strong infectious or hub model
isHub = true
# hub constant
hubconstant = 6^0.5

#### sets to store different compartments#####
# Susceptible set
S = Set()
# exposed set
E = Set()
# infectious but not ICU set
I = Set()
# infectious and going to go to ICU set
L = Set()
# ICU set
ICU = Set()
# removed set
R = Set()
# current ICU capacity of NYC according to NYPost
# 2.7 beds per 10000 according to NYPost
capacity = (2.7/10000) * N
### arrays for plotting purposes
plotS = zeros(t+1)
# exposed
plotE = zeros(t+1)
# infected
plotI = zeros(t+1)
# lag time people
plotL = zeros(t+1)
# ICU
plotICU = zeros(t+1)
# Removed people
plotR = zeros(t+1)
# infectious people array
plotInfectious = zeros(t+1)
# array for time
tvec = zeros(t+1)
for i in 2:t+1
    tvec[i] = i-1
end
# find the distance between two people
function distance(p1, p2)
    x1 = p1[1]
    x2 = p2[1]
    y1 = p1[2]
    y2 = p2[2]
    return (((x1-x2)^2) + ((y1-y2)^2))^0.5
end
# generate a boolean value to represent if an event takes place
function randomEvent(p :: Float64)
    u = rand(Float64)
    if(p == 0)
        return false
    elseif (0<u<=p)
        return true
    else
        return false
    end
end
# generate the probability

function initializePop()
# determine if the first infected person is a superspreader
    ss1 = randomEvent(λ)
    first = [sizeL/2.0, 0.0, ss1]
    push!(I, first)
    # generate the rest of the population
    for i in 1:N-1
        ss = randomEvent(λ)
        x = rand(Float64) * sizeL
        y = rand(Float64) * sizeL
        p = [x,y,ss]
        push!(S, p)
    end
    plotS[1] = length(S)
    plotI[1] = length(I)
    plotInfectious[1] = plotI[1] + plotL[1]
end

function infectProb4(r :: Float64, ss :: Bool)
    # if not a super spreader
    if ss == false
        if r > rstart
            return 0.0
        else
            w = w0 * (1.0 - (r/rstart)^α)
            return w
        end
    # is a super spreader
    else
        # if it is is the hub model
        if isHub == true
            rnew = hubconstant * rstart
            if r > rnew
                return 0.0
            else
                return w0 * (1.0-(r/rnew)^α)
            end
        # if it is the strong infectious model
        else
            if r > rstart
                return 0
            else
                return w0
            end
        end
    end
end
# runs infectious infecting S to E
function infect(tempInf, tempSus)

    for p in tempInf
        for q in tempSus
            d = distance(p,q)
            w = infectProb4(d, convert(Bool,p[3]))
            event = randomEvent(w)
            if event == true
                delete!(S, q)
                push!(E, q)
            end
        end
    end
end
# run Exposed transition
function E_IL(tempE)
    for p in tempE
        # check if exposed becomes infectious
        event = randomEvent(ζ)
        if event == true
            delete!(E, p)
            # check if someone goes L(event2 = true) or I
            event2 = randomEvent(ρ)
            if event2 == true
                push!(L, p)
            else
                push!(I, p)
            end
        end
    end
end
# run simulation of the recovery
function recovery(tempI, tempICU)
    for p in tempI
        event = randomEvent(γ)
        if event == true
            delete!(I, p)
            push!(R, p)
        end
    end
    for p in tempICU
        event = randomEvent(κ)
        if event == true
            delete!(ICU, p)
            push!(R, p)
        end
    end
end
# L to ICU
function LtoICU(tempL)
    for p in tempL
        event = randomEvent(η)
        if event == true
            delete!(L, p)
            push!(ICU, p)
        end
    end
end

# run state change
function simulate4(pos)
    infectious = union!(I, L)
    tempS = S
    tempE = E
    tempI = I
    tempL = L
    tempICU = ICU
    # run state change for infectious people
    infect(infectious, tempS)
    # run state change from E to I or L
    E_IL(tempE)
    # runs recovered state change
    recovery(tempI, tempICU)
    # runs L to ICU
    LtoICU(tempL)
    # put the lengths of sets in plottable matrices
    plotS[pos] = length(S)
    plotE[pos] = length(E)
    plotI[pos] = length(I)
    plotL[pos] = length(L)
    plotICU[pos] = length(ICU)
    plotInfectious[pos] = plotI[pos] + plotL[pos]
    plotR[pos] = length(R)
end
initializePop()
# run multiple state changes
function run2()
    for i in 2:t+1
        simulate4(i)
    end
end
function visualize()

    #plot!(p, tvec, plotInfectious, label = "Infectious Individuals")
    plot!(tvec, plotICU, label = "run3" )
    xlabel!("Number of Days")
    ylabel!("Number of Hospitalizations")
    title!("ICU Capacity")
end

run2()
visualize()
savefig("N=100000 ICU Runs")
