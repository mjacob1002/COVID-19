using DataFrames, CSV
#=
S_df = DataFrame(); I_df = DataFrame(); R_df = DataFrame()
E_df = DataFrame(); L_df = DataFrame(); ICU_df = DataFrame()
=#


function load_data(out, name)
    df = DataFrame()
    nms = [Symbol("t$i") for i in 1:length(out)]
    for j in 1:length(out)
        df[nms[j]] = out[j]
    end
    CSV.write("/Users/mathewjacob/Desktop/RxCOVea/COVID-19/SpatialModel/MonteCarlo/Results/GoodResults/ShutDown/"*name*".csv", df)
    return df
end

S_df = load_data(S_out, "S_df")
I_df = load_data(I_out, "I_df")
R_df = load_data(R_out, "R_df")
E_df = load_data(E_out, "E_df")
L_df = load_data(L_out, "L_df")
ICU_df = load_data(ICU_out, "ICU_df")

nms = [Symbol("t$i") for i in 1:num_MC_simulations]
total_infected = DataFrame()
for i in 1:num_MC_simulations
    curr = Any[]
    for j in 1:480
        I_curr = I_df[nms[i]][j]
        L_curr = L_df[nms[i]][j]
        ICU_curr = ICU_df[nms[i]][j]

        curr_total_infected = I_curr + L_curr + ICU_curr
        push!(curr, curr_total_infected)
    end
    total_infected[nms[i]] = curr
end

CSV.write("/Users/mathewjacob/Desktop/RxCOVea/COVID-19/SpatialModel/MonteCarlo/Results/GoodResults/ShutDown/total_infected.csv", total_infected)

#DF-->CSV
