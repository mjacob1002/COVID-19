using Random, Distributions
using CoordinateTransformations, StaticArrays

include("./0_ICU_circular_mobility_functions.jl")
#Parameters
num_MC_simulations = 10

N = 10000#150-900
sampleSize = 200
loc_lmean = 450.0 #no unit
loc_lstd = 40
w0 = 0.9 #infection probability when dist =0.0
λmean = .15 #ratio of superspreaders
γmean = 1/10 #rate of recovery and death
ζmean = .33 #E-->Infection
ρmean = .04 #E-->L # probability of, once infected, going to L compartment
ηmean = .3 #L-->ICU
κmean = 1/14 #ICU-->R
tf = 40 #time unit unspecified (can be day or hour or even second) #40 # the number of full periods
model = "hub" #"hub" or "strong infectiousness"
time_steps = 12.0 #determines periodical movement e.g 1 period every 12 days
c_from_p = CartesianFromPolar()
p_from_c = PolarFromCartesian()
const global delta = (2 * π)/time_steps
hospital_x = loc_lmean * 1.5
hospital_y = loc_lmean * 0.5

S_out = Any[]; I_out = Any[]; R_out = Any[]
E_out = Any[]; L_out = Any[]; ICU_out = Any[]

@time begin
for a in 1:num_MC_simulations
    # initialize the variables from distribution
    λ = rand(Normal(λmean, stdDev(λmean,sampleSize)))
    # nature.com/articles/s41598 "Efficacy of Masks and Face Coverings"
    r0 = rand(Normal(3, 0.5)) #OG: Normal(3,0.5)
    γ = rand(Normal(γmean, stdDev(γmean, sampleSize)))
    ζ = rand(Normal(ζmean, stdDev(ζmean,sampleSize)))
    ρ = rand(Normal(ρmean, stdDev(ρmean, sampleSize)))
    η = rand(Normal(ηmean, stdDev(ηmean, sampleSize)))
    κ = rand(Normal(κmean,stdDev(κmean,sampleSize)))
    loc_l = rand(Normal(loc_lmean, loc_lstd))
    μ = 45/4 # OG 45
    σ = 10/4 # OG 10
    @time begin
    println("run# "*string(a))
    println("ρ:", ρ, "ζ: ", ζ, " r0: ", r0)
    pop = initialize_pop(loc_l, λ, N)

    pop = generate_random_vector(pop, N)
    pop = generate_random_radii(pop, μ, σ, N)
    pop = genetate_circular_path(pop)

    S, E, L, ICU, I, R = parse_pop(pop)

    #initialize output vectors -->output as dataframes later
    S_vals = Int64[]; I_vals = Int64[]; R_vals = Int64[]
    E_vals = Int64[]; L_vals = Int64[]; ICU_vals = Int64[]

    for b in 1:tf

        S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals, S, E, L, ICU, I, R = move(S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals, S, E, L, ICU, I, R, N, r0, ρ,ζ,η,κ,γ) #e.g. 40 days
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
