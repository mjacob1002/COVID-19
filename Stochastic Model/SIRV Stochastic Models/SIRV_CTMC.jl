using Random
using Plots

#parameters
S0 = 995
I0 = 5
R0 = 0
V0 = 0
beta = 1/100
gamma = 1/1000
p = 1/10 # vaccination rate
N = S0 + I0 + R0 + V0
SIV = [S0 I0 V0] #state 0

#time and state arrays for plotting
t_vec = Float64[]
S_vals = Int[]
I_vals = Int[]
V_vals = Int[]
push!(t_vec, 0)
push!(S_vals, S0)
push!(I_vals, I0)
push!(V_vals, V0)

#simulation function based on CTMC
#based on Allen, 2017
function gillespie_simulation(SIV)
    s = SIV[1] #get current s and i
    i = SIV[2]
    v = SIV[3]
    lambda = (beta * s * i / N) + (gamma * i) + (p * s / N)

    #u1 - inter-event time, U[0,1]
    u1 = rand(Float64)
    tao = - log(u1)/lambda #(6) in Allen, 2017
    dt = tao
    # u2 - event, U[0,1]
    u2 = rand(Float64)
    p1 = (beta * s * i / N)/lambda
    p2 = (gamma * i)/lambda
    p3 = (p * s / N)/lambda

    #println(u2, p1, p2)
    #FIXME: number of states???
    if u2 >=0 && u2 < p1 #event 1: infection (-1, +1, 0)
        SIV[1] = SIV[1] - 1
        SIV[2] = SIV[2] + 1
    elseif u2 >= p1 && u2 < p1+p2 #event 2: recovery (0, -1, 0)
        SIV[2] = SIV[2] - 1
    elseif u2 >= p1+p2 && u2 < p1+p2+p3 #event 3: vaccination (-1, 0, +1)
        SIV[1] = SIV[1] - 1
        SIV[3] = SIV[3] + 1
    elseif u2 >= p1+p2+p3 #no change?

    else #otherwise?

    end
    return dt, SIV
end

#main loop: simulation
while SIV[2] != 0
    #note scope issue: modifying global variable
    global dt, SIV = gillespie_simulation(SIV)
    push!(t_vec, dt)
    push!(S_vals, SIV[1])
    push!(I_vals, SIV[2])
    push!(V_vals, SIV[3])
end

#dt-->continuous t_vec for plotting
time_elapsed = 0
for (index, value) in enumerate(t_vec)
    global time_elapsed = time_elapsed + value
    global t_vec[index] = time_elapsed
end

#add Removed
R_vals = Int[]
for (index, value) in enumerate(S_vals)
    r = N - S_vals[index] - I_vals[index] - V_vals[index]
    push!(R_vals, r)
end
print(S_vals)
#plot
p = plot(t_vec, S_vals,
    title = "SIRV CTMC simulation using Gillespie algorithm",
    label = "Susceptible")
plot!(p, t_vec, I_vals,
    label = "Infected")
plot!(p, t_vec, V_vals,
    label = "Vaccinated")
plot!(p, t_vec, R_vals,
    label = "Removed")
