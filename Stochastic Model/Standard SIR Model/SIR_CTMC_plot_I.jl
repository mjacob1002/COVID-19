#SIR: CTMC vs deterministic

using Random
using Plots

#parameters
S0 = 98 #748, 98 from Allen et al, 2017
I0 = 2
R0 = 0
beta = 0.3 #dynamic-->hammer and dance
gamma = 0.15
N = S0 + I0 + R0
SI = [S0 I0]
time_elapsed = 0.0

#time and state arrays for plotting
t_vec_1 = Float64[] #stores time stamps
t_vec_2 = Float64[]
t_vec_3 = Float64[]
t_vec_4 = Float64[]
I_vals_1 = Int[]
I_vals_2 = Int[]
I_vals_3 = Int[]
I_vals_4 = Int[]
#initialize
push!(t_vec_1, 0.0)
push!(t_vec_2, 0.0)
push!(t_vec_3, 0.0)
push!(t_vec_4, 0.0)
push!(I_vals_1, I0)
push!(I_vals_2, I0)
push!(I_vals_3, I0)
push!(I_vals_4, I0)

#simulation function using CTMC
#based on Allen, 2017
function gillespie_simulation(run, SI)
    s = SI[1] #get current s and i
    i = SI[2]
    lambda = (beta * s * i / N) + gamma * i

    #u1 - inter-event time, U[0,1]
    u1 = rand(Float64)
    tao = - log(u1)/lambda #(6) in Allen, 2017
    dt = tao
    global time_elapsed = time_elapsed + dt
    if run == 1
        push!(t_vec_1, time_elapsed)
    elseif run == 2
        push!(t_vec_2, time_elapsed)
    elseif run == 3
        push!(t_vec_3, time_elapsed)
    elseif run == 4
        push!(t_vec_4, time_elapsed)
    end


    # u2 - event, U[0,1]
    u2 = rand(Float64)
    p1 = (beta * s * i / N)/lambda #probability of infection
    p2 = (gamma * i)/lambda #probability of recovery

    #println(u2, p1, p2)
    if u2 >=0 && u2 < p1 #event 1: infection
        SI[1] = SI[1] - 1
        SI[2] = SI[2] + 1
    elseif u2 >= p1 && u2 < p1+p2 #event 2: recovery
        SI[2] = SI[2] - 1
    elseif u2 >= p1+p2 #no change?

    else #otherwise?

    end
    return SI
end


#main loop: simulation
#4 simulations
for i = 1:4
    global SI = [S0 I0]#state 0
    global time_elapsed = 0.0
    while SI[2] != 0
        global SI = gillespie_simulation(i, SI)
        if i == 1
            push!(I_vals_1, SI[2])
        elseif i == 2
            push!(I_vals_2, SI[2])
        elseif i == 3
            push!(I_vals_3, SI[2])
        elseif i == 4
            push!(I_vals_4, SI[2])
        end
    end
end


########################################################
#deterministic SIR model
function updateSIR(SI)
    s = SI[1]
    i = SI[2]
    s_new = s - beta * i * s * dt / N
    i_new = i + beta * i * s * dt / N - gamma * i * dt
    return [s_new, i_new]
end

t_final = 100. # days
dt = 0.5
n_steps = round(Int64, t_final/dt) # round function to ensure integer n-steps
t_vec = Float64[] # initialise array for time
push!(t_vec, 0.0)
I_val = Float64[]
push!(I_val, I0)
SI = [S0 I0]
for j = 1:n_steps
    global SI = updateSIR(SI)
    push!(I_val, SI[2])
    push!(t_vec, t_vec[j]+dt)
end
#println(t_vec)

#plot
p = plot(t_vec_1, I_vals_1,
    title = "SIR: CTMC vs deterministic",
    label = "run1")
plot!(p, t_vec_2, I_vals_2,
    label = "run2")
plot!(p, t_vec_3, I_vals_3,
    label = "run3")
plot!(p, t_vec_4, I_vals_4,
    label = "run4")
plot!(p, t_vec, I_val,
    label = "deterministic")
xlabel!("days")
ylabel!("I(t)")
