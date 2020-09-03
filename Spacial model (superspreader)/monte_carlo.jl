using Random
using Plots
using DataFrames, CSV

include("ss_functions.jl")

#---------Parameters--------
N = 200 #150-900
r0 = 2.0
L = 10 * r0
w0 = 1.0
γ = 0.1 #recovery or death
λ = C_NULL #ss density#0.0-->1.0
λ_vals = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
model = "hub"#"hub" or "strong infectiousness"
#---------Simulation--------------

results = Any[]
I = Any[]
S_vals = Int[]
I_vals = Int[]
R_vals = Int[]
Susceptible = Set()
Infected = Set()
Removed = Set()
for i in 1:6
    global λ = λ_vals[i]
    #println(λ)
    for j in 1:100 #100 monte carlo simulations
        initialize_run(N, L) # a set of all ppl
        push!(S_vals, length(Susceptible))
        push!(I_vals, length(Infected))
        push!(R_vals, length(Removed))
        for k in 1:40
            simulate()
            push!(S_vals, length(Susceptible))
            push!(I_vals, length(Infected))
            push!(R_vals, length(Removed))
        end
        push!(I, I_vals)
        global S_vals = Int[]
        global I_vals = Int[]
        global R_vals = Int[]
        global Susceptible = Set()
        global Infected = Set()
        global Removed = Set()
    end
    push!(results, I)
    global I = Any[]
end

#simulations (I(t)) stored in results -->plotting/ analysis
