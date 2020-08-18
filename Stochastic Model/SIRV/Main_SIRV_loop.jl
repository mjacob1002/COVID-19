# Title: SIR model example, looking at dynamics over days
#
# Name: Jantine Broek
# Date: March 20
#vaccination added by Eva Xueyao Guo on June 20 2020
##############################################################

# working directory
#cd("//Users/jantinebroek/Documents/03_prog_comp/02_collab/COVID")

# Loads packages
using Plots
pyplot()    # Plots package will use Pyplot (matplotlib needs to be installed).

# Functions
include("updateSIRV.jl")

# change parameters here
λ = 1/1000
γ = 1/100
p = 1/100

# input
s0 = 98.
i0 = 2.
r0 = 0.
v0 = 0.

# time steps
t_final = 610 # days
dt = 0.5
n_steps = round(Int64, t_final/dt) # round function to ensure integer n-steps

# pre-allocate
t_vec = Array{Float64}(undef, n_steps + 1) # initialise array for time
result_vals = Array{Float64}(undef, n_steps + 1, 4) # initialise array for results

# assign
result_vals[1,:] = [s0, i0, r0, v0]
t_vec[1] = 0

for j = 1:n_steps
    # call function
    result_vals[j + 1, :] = updateSIR(result_vals[j, :])
    t_vec[j + 1] = t_vec[j] + dt
end

# plot
plot(t_vec, result_vals,
    title = "SIRV results",
    xlabel = "days",
    ylabel = "num ppl",
    label = ["Susceptible" "Infected" "Removed" "Vaccinated"]
    )

#savefig("/Users/xueyaoguo/Desktop/SIRV_result.png")
