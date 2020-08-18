# Title: SIR model example, looking at dynamics over days
#
# Name: Jantine Broek
# Date: March 20
# vaccination added by Eva Xueyao Guo on June 20 2020
##############################################################

# working directory
#cd("//Users/jantinebroek/Documents/03_prog_comp/02_collab/COVID")

# Loads packages

# Functions
include("updateSIRV.jl")

#parameters
dt = 0.5
λ = 1/2000
γ = 1/1000
p = 1/10

# input
s, i, r, v= 1000., 10., 20. ,0.  # multiple assignments
pop_vec = [s i r v]

# call function
(S_up, I_up, R_up, V_up) = updateSIR(pop_vec)
