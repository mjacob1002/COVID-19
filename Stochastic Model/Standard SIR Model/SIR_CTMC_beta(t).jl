#TODO: beta(t) to explore hammer and dance

using Random
using Plots

#parameters
S0 = 900
I0 = 100
R0 = 0
beta0 = 0.3 #initial beta value
beta = 0.3 #updated in every step
beta_scaling_factor = 1 #to scale beta, smaller scaling factor-->smaller effect
time_elapsed = 0.0 #updated in every step
gamma = 0.15
N = S0 + I0 + R0
SI = [S0 I0] #state 0

#time and state arrays for plotting
t_vec = Float64[]
S_vals = Int[]
I_vals = Int[]
push!(t_vec, 0.0)
push!(S_vals, S0)
push!(I_vals, I0)

function update_beta_t()
    #continuous time
    #simple oscilation examlpe for beta(t)
    #can also implement as:
        # a different function
        # real beta value measured from data
    global beta = sin(time_elapsed) * beta_scaling_factor
end


#simulation function based on CTMC
#from Allen, 2017
function gillespie_simulation(SI)
    s = SI[1] #get current s and i
    i = SI[2]
    update_beta_t()
    lambda = (beta0 * s * i / N) + gamma * i

    #u1 - inter-event time, U[0,1]
    u1 = rand(Float64)
    tao = - log(u1)/lambda #(6) in Allen, 2017
    dt = tao
    global time_elapsed = time_elapsed + dt #update global time elapsed
    #println(dt)
    push!(t_vec, time_elapsed)

    # u2 - event, U[0,1]
    u2 = rand(Float64)
    p1 = (beta * s * i / N)/lambda #probability of infection
#    p2 = (gamma * i)/lambda #probability of recovery
    #Note that different implementation here will affect the result
    p2 = 1 - p1
    #p3 = 1 - p1 - p2

    #println(u2, p1, p2)
    if u2 >=0 && u2 < p1 #event 1: infection
        SI[1] = SI[1] - 1
        SI[2] = SI[2] + 1
    elseif u2 >= p1 && u2 < p1+p2 #event 2: recovery
        SI[2] = SI[2] - 1
#    elseif u2 >= p1+p2 && u2 < p1+p2+p3 #no change

    else #otherwise

    end
    return SI
end

#main loop: simulation
while SI[2] != 0
    #note scope issue: modifying global variable
    global SI = gillespie_simulation(SI)
    push!(S_vals, SI[1])
    push!(I_vals, SI[2])
end

#add Removed for plotting
R_vals = Int[]
for (index, value) in enumerate(S_vals)
    r = N - S_vals[index] - I_vals[index]
    push!(R_vals, r)
end

#plot
p = plot(t_vec, S_vals,
    title = "SIR CTMC with β(t)=sin(t)+β0",
    label = "Susceptible")
plot!(p, t_vec, I_vals,
    label = "Infected")
plot!(p, t_vec, R_vals,
    label = "Removed")
xlabel!("days")
ylabel!("num ppl")
