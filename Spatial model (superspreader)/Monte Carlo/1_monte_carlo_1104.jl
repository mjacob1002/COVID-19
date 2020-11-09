using Random, Distributions
using CoordinateTransformations, StaticArrays

include("./ICU_circular_mobility_functions.jl")

#Parameters
num_MC_simulations = 100

N = 1000 #150-900
loc_l = 100.0 #no unit
r0 = loc_l/5.0
w0 = 1.0 #infection probability when dist =0.0
λ = 1.0 #ratio of superspreaders
γ = 0.05 #rate of recovery and death
ζ = .5 #E-->I
ρ = .03 #E-->L # probability of, once infected, going to L compartment
η = .8 #L-->ICU
κ = .3 #ICU-->R
tf = 40 #time unit unspecified (can be day or hour or even second) #40
model = "hub" #"hub" or "strong infectiousness"
time_steps = 12.0 #determines periodical movement
μ = loc_l
σ = loc_l
c_from_p = CartesianFromPolar()
p_from_c = PolarFromCartesian()
const global delta = (2 * π)/time_steps
hospital_x = loc_l * 1.5
hospital_y = loc_l * 0.5

S_out = Any[]; I_out = Any[]; R_out = Any[]
E_out = Any[]; L_out = Any[]; ICU_out = Any[]

@time begin
for a in 1:num_MC_simulations
    @time begin
    println("run# "*string(a))
    pop = initialize_pop(loc_l, λ, N)

    pop = generate_random_vector(pop, N)
    pop = generate_random_radii(pop, μ, σ, N)
    pop = genetate_circular_path(pop)

    S, E, L, ICU, I, R = parse_pop(pop)

    #initialize output vectors -->output as dataframes later
    S_vals = Int64[]; I_vals = Int64[]; R_vals = Int64[]
    E_vals = Int64[]; L_vals = Int64[]; ICU_vals = Int64[]

    for b in 1:tf

        S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals, S, E, L, ICU, I, R = move(S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals, S, E, L, ICU, I, R, N) #e.g. 40 days
        #=
        #after each day/time unit: output num ppl in each group
        push!(S_vals, length(S)); push!(I_vals, length(I)); push!(R_vals, length(R))
        push!(E_vals, length(E)); push!(L_vals, length(L)); push!(ICU_vals, length(ICU))
        =#
    end
    end

    #after each simulation round: output time series simulation results
    push!(S_out, S_vals); push!(I_out, I_vals); push!(R_out, R_vals)
    push!(E_out, E_vals); push!(L_out, L_vals); push!(ICU_out, ICU_vals)
end
println("DONE")

end
#save to a dataframe
