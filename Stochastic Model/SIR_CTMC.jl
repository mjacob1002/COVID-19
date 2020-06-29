using Random
using Plots

#parameters
S0 = 995
I0 = 5
R0 = 0
beta = 1/100
gamma = 1/1000
N = S0 + I0 + R0
SI = [S0 I0] #state 0

#time and state arrays for plotting
t_vec = Float64[]
S_vals = Int[]
I_vals = Int[]
push!(t_vec, 0)
push!(S_vals, S0)
push!(I_vals, I0)

#simulation function using CTMC
#based on Allen, 2017
function gillespie_simulation(SI)
    s = SI[1] #get current s and i
    i = SI[2]
    lambda = (beta * s * i / N) + gamma * i

    #u1 - inter-event time, U[0,1]
    u1 = rand(Float64)
    tao = - log(u1)/lambda #(6) in Allen, 2017
    dt = tao
    # u2 - event, U[0,1]
    u2 = rand(Float64)
    p1 = (beta * s * i / N)/lambda
    p2 = (gamma * i)/lambda

    #println(u2, p1, p2)
    #FIXME: number of states???
    if u2 >=0 && u2 < p1 #event 1: infection
        SI[1] = SI[1] - 1
        SI[2] = SI[2] + 1
    elseif u2 >= p1 && u2 < p1+p2 #event 2: recovery
        SI[2] = SI[2] - 1
    elseif u2 >= p1+p2 #no change?

    else #otherwise?

    end
    return dt, SI
end

#main loop: simulation
while SI[2] != 0
    #note scope issue: modifying global variable
    global dt, SI = gillespie_simulation(SI)
    push!(t_vec, dt)
    push!(S_vals, SI[1])
    push!(I_vals, SI[2])
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
    r = N - S_vals[index] - I_vals[index]
    push!(R_vals, r)
end

#plot
p = plot(t_vec, S_vals,
    title = "SIR CTMC simulation using Gillespie algorithm",
    label = "Susceptible")
plot!(p, t_vec, I_vals,
    label = "Infected")
plot!(p, t_vec, R_vals,
    label = "Removed")
