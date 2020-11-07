using DataFrames, CSV
using Statistics
using Plots

S = CSV.read(joinpath(Pkg.dir("DataFrames"), "/Users/xueyaoguo/Desktop/Monte_Carlo/lamda_1.0/S_df.csv"));
I = CSV.read(joinpath(Pkg.dir("DataFrames"), "/Users/xueyaoguo/Desktop/Monte_Carlo/lamda_1.0/I_df.csv"));
R = CSV.read(joinpath(Pkg.dir("DataFrames"), "/Users/xueyaoguo/Desktop/Monte_Carlo/lamda_1.0/R_df.csv"));
E = CSV.read(joinpath(Pkg.dir("DataFrames"), "/Users/xueyaoguo/Desktop/Monte_Carlo/lamda_1.0/E_df.csv"));
L = CSV.read(joinpath(Pkg.dir("DataFrames"), "/Users/xueyaoguo/Desktop/Monte_Carlo/lamda_1.0/L_df.csv"));
ICU = CSV.read(joinpath(Pkg.dir("DataFrames"), "/Users/xueyaoguo/Desktop/Monte_Carlo/lamda_1.0/ICU_df.csv"));
total_infected = CSV.read(joinpath(Pkg.dir("DataFrames"), "/Users/xueyaoguo/Desktop/Monte_Carlo/lamda_1.0/total_infected.csv"));

nms = [Symbol("t$i") for i in 1:num_MC_simulations]
total_IR = DataFrame()
for i in 1:100
    curr = Any[]
    for j in 1:480
        I_curr = I[nms[i]][j]
        L_curr = L[nms[i]][j]
        ICU_curr = ICU[nms[i]][j]
        R_curr = R[nms[i]][j]

        curr_total_IR = I_curr + L_curr + ICU_curr + R_curr
        push!(curr, curr_total_IR)
    end
    total_IR[nms[i]] = curr
end

#CSV.write("/Users/xueyaoguo/Desktop/Monte_Carlo/lamda_1.0/total_IR.csv", total_IR)

#DATAFRAME total_IR

N = 1000
num_MC_simulations = 100
nms = [Symbol("t$i") for i in 1:num_MC_simulations]
IR_max = Any[]
#find the date for 10%, 20%, 30%, ...
percent_time = Any[]
for i in 1:num_MC_simulations
    curr_time = Int32[]
    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.10
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.20
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.30
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.40
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.50
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.60
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.70
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.80
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 0.90
            push!(curr_time, j)
            break
        end
    end

    for j in 1:480 #time
        if total_IR[nms[i]][j]/N >= 1.00
            push!(curr_time, j)
            break
        end
    end
    push!(percent_time, curr_time)
end

#take average
time_mean = Float64[]
for i in 1:10
    curr_mean = 0.0
    for j in 1:num_MC_simulations
        curr_mean = curr_mean + percent_time[j][i]
    end
    curr_mean = curr_mean/Float64(num_MC_simulations)
    push!(time_mean, curr_mean)
end

println(time_mean)

#and do it for all lamda values
