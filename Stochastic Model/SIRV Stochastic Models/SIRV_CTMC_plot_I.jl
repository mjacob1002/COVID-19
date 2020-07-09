#SIRV: CTMC vs deterministic

using Random
using Plots

#parameters
S0 = 98
I0 = 2
R0 = 0
V0 = 0
beta = 0.3 #dynamic-->hammer and dance
gamma = 0.15
p = 1/10 # vaccination rate
N = S0 + I0 + R0 + V0
SIV = [S0 I0 V0] #state 0
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
function gillespie_simulation(run, SIV)
    s = SIV[1] #get current s and i
    i = SIV[2]
    v = SIV[3]
    lambda = (beta * s * i / N) + (gamma * i) + (p * s / N)

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
    p3 = (p * s / N)/lambda #probability of vaccination

    #println(u2, p1, p2)
    if u2 >=0 && u2 < p1 #event 1: infection (-1, +1, 0)
        SIV[1] = SIV[1] - 1
        SIV[2] = SIV[2] + 1
    elseif u2 >= p1 && u2 < p1+p2 #event 2: recovery (0, -1, 0)
        SIV[2] = SIV[2] - 1
    elseif u2 >= p1+p2 && u2 < p1+p2+p3 #event 3: vaccination (-1, 0, +1)
        SIV[1] = SIV[1] - 1
        SIV[3] = SIV[3] + 1
    elseif u2 >= p1+p2+p3 #no change

    else #otherwise

    end
    return SIV
end


#main loop: simulation
#4 simulations

for i = 1:4
    global SIV = [S0 I0 V0]#state 0
    global time_elapsed = 0.0
    while SIV[2] != 0
        global SIV = gillespie_simulation(i, SIV)
        if i == 1
            push!(I_vals_1, SIV[2])
        elseif i == 2
            push!(I_vals_2, SIV[2])
        elseif i == 3
            push!(I_vals_3, SIV[2])
        elseif i == 4
            push!(I_vals_4, SIV[2])
        end
    end
end


########################################################
#deterministic SIR model
function updateSIR(SIV)
    s = SIV[1]
    i = SIV[2]
    v = SIV[3]
    s_new = s - beta * i * s * dt / N - p * s * dt / N
    i_new = i + beta * i * s * dt / N - gamma * i * dt
    v_new = v + p * s * dt / N
    return [s_new, i_new, v_new]
end

t_final = 100. # days
dt = 0.5
n_steps = round(Int64, t_final/dt) # round function to ensure integer n-steps
t_vec = Float64[] # initialise array for time
push!(t_vec, 0.0)
I_val = Float64[]
push!(I_val, I0)
SIV = [S0 I0 V0]
for j = 1:n_steps
    global SIV = updateSIR(SIV)
    push!(I_val, SIV[2])
    #t_vec[j + 1] = t_vec[j] + dt
    push!(t_vec, t_vec[j]+dt)
end
println(t_vec)
#plot
p = plot(t_vec_1, I_vals_1,
    title = "SIRV: CTMC vs deterministic",
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
