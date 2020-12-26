using DataFrames, CSV
using Statistics
using Plots

ICU = CSV.read(joinpath(Pkg.dir("DataFrames"), "/Users/mathewjacob/Desktop/RxCOVea/COVID-19/SpatialModel/MonteCarlo/Results/GoodResults/ShutDown/ICU_df.csv"));
N = 10000
num_MC_simulations = 10
#x = [0:100:480]
fig = plot()
nms = [Symbol("t$i") for i in 1:num_MC_simulations]

for i in 1:num_MC_simulations
    plot!(fig, ICU[nms[i]][1:480], label="run "*string(i))
end

xlabel!(fig, "days")
ylabel!(fig, "num ICU patients")
title!(fig, "ICU Patients When Some Businesses Shut Down")
display(fig)
savefig("/Users/mathewjacob/Desktop/RxCOVea/COVID-19/SpatialModel/MonteCarlo/Results/GoodResults/ShutDown/shutdown.png")
