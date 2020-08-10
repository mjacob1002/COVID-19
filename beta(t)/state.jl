using Plots
cd("/Users/evaxueyaoguo/Desktop/results/")
#make sure you have a directory called "New York"
#i.e. make /Users/evaxueyaoguo/Desktop/results/New York
#and /Users/evaxueyaoguo/Desktop/results/New York/figs
#before running code
state_name = "New York" #TODO: change this to any state name
#mkdir(state_name)

state = Any[]
n = length(group_by_county)
for i = 1:n
    if group_by_county[i].state[1]== state_name
        curr_df = group_by_county[i]
        push!(state, curr_df)
    end
end

#US cumulative recovery rate calculated from JHU global time series data
#see SIR_US.xls for recovery rate calculations
#2020-01-21 to 2020-08-03
#cumulative_recovered, 0 added for 2020-01-21
cumulative_recovery_rate = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.272727273,0.272727273,0.25,0.25,0.230769231,0.230769231,0.230769231,0.230769231,0.230769231,0.230769231,0.230769231,0.230769231,0.333333333,0.333333333,0.333333333,0.333333333,0.4,0.4,0.375,0.4375,0.291666667,0.233333333,0.132075472,0.095890411,0.067307692,0.040229885,0.031531532,0.020771513,0.015521064,0.013487476,0.011251758,0.007213706,0.00768738,0.005563282,0.004181185,0.004043127,0.003899083,0.002768279,0.011775261,0.008547009,0.007546589,0.006815102,0.005272356,0.004059293,0.006431106,0.005465143,0.00809837,0.008496617,0.008781918,0.018873269,0.03468812,0.037218372,0.039560234,0.036797351,0.035100724,0.047310453,0.051686598,0.053322985,0.054682004,0.05482841,0.054710814,0.057817863,0.059226962,0.059275285,0.074735353,0.078444286,0.0816585,0.081735079,0.083397911,0.088780933,0.092992233,0.092289809,0.092922683,0.092385918,0.092016333,0.109128899,0.106684637,0.110465681,0.112437852,0.114164397,0.115738832,0.143517979,0.148184589,0.154382302,0.155088063,0.158079734,0.157076517,0.153956906,0.154617574,0.154427291,0.161706434,0.162035524,0.172017396,0.167462942,0.174363052,0.17307672,0.173045085,0.18205746,0.182344896,0.186928387,0.188486013,0.188763543,0.188334326,0.217657257,0.221554319,0.222090742,0.22700206,0.227865517,0.229294648,0.231173928,0.23162431,0.234070213,0.247208086,0.252263307,0.252462482,0.25803572,0.258160637,0.258261564,0.259959961,0.260427285,0.264311329,0.265090065,0.26665787,0.266988065,0.267149702,0.268303076,0.268251108,0.272623894,0.272954361,0.273745545,0.273431278,0.272977923,0.27377836,0.272654044,0.276865964,0.275846851,0.275417159,0.273938932,0.271851801,0.270612714,0.268706096,0.272208944,0.273337571,0.271616781,0.285177252,0.282755608,0.314765625,0.313636842,0.31475605,0.31256521,0.312129608,0.310817121,0.308620894,0.306715651,0.304491274,0.306745196,0.305719183,0.307491322,0.304976823,0.303533582,0.302504733,0.299772875,0.30252535,0.303142867,0.304990452,0.305354094,0.30677556,0.306155344,0.306539113,0.309020946,0.311133817,0.313853772,0.314605179,0.315245072,0.316394918,0.31463221,0.32108479]
crr_dates = collect(Date(2020,1,21):Day(1):Date(2020,8,3))
rate_df = DataFrame()
rate_df.date = crr_dates
rate_df.cumulative_recovery_rate = cumulative_recovery_rate

#get number of counties
num_counties = length(state)

for i = 1:num_counties #for every county in this state

    #get name of current county
    curr_county = state[i].county[1]
    state_dates = state[i].date
    println(curr_county)

    state_df = DataFrame()
    state_df.date = state_dates
    state_df.cumulative_cases = state[i].cases
    state_df.cumulative_deaths = state[i].deaths

    curr_df = innerjoin(rate_df, state_df, on= :date)

    cumulative_cases = curr_df.cumulative_cases
    #calculate new cases
    new_cases = Int64[]
    push!(new_cases, cumulative_cases[1])
    n = length(cumulative_cases)
    #print(n)
    for i = 2:n
        curr = cumulative_cases[i]- cumulative_cases[i-1]
        push!(new_cases, curr)
    end

    cumulative_recovery = Int64[]
    push!(cumulative_recovery, round(curr_df.cumulative_recovery_rate[1]*cumulative_cases[1]))
    for i = 2:n
        push!(cumulative_recovery, round(curr_df.cumulative_recovery_rate[i]*cumulative_cases[i]))
    end

    new_recovery = Int64[]
    push!(new_recovery, cumulative_recovery[1])
    for i = 2:n
        curr = cumulative_recovery[i]-cumulative_recovery[i-1]
        if curr < 0 #correct negative values (approximation only)
            curr = 0
        end
        push!(new_recovery, curr)
    end

    cumulative_death = curr_df.cumulative_deaths
    new_death = Int64[]
    push!(new_death, cumulative_death[1])
    for i = 2:n
        curr = cumulative_death[i]-cumulative_death[i-1]
        if curr < 0 #correct negative values (approximation only)
            curr = 0
        end
        push!(new_death, curr)
    end

    #R(t), removed
    R_vals = Int64[] #R(t) + recovered + death
    push!(R_vals, new_recovery[1] + new_death[1])
    for i = 2:n
        push!(R_vals, R_vals[i-1] + new_recovery[i] + new_death[i])
    end
    #print(R_vals)

    #I(t), infected
    I_vals = Int64[] #I(t) + new - recovered + death
    push!(I_vals, new_cases[1] - new_recovery[1] - new_death[1])
    for i = 2:n
        push!(I_vals, I_vals[i-1] + new_cases[i] - new_recovery[i] - new_death[i])
    end
    #print(I_vals)

    #β(t) and γ(t) from Chen et al.
    beta_vals = Float64[]
    gamma_vals = Float64[]
    R_0 = Float64[]
    for i = 1:n-1
        #eq. 11 and 12 from Chen et al.
        beta = ((I_vals[i+1]-I_vals[i])+(R_vals[i+1]-R_vals[i]))/(I_vals[i])
        gamma = (R_vals[i+1]-R_vals[i])/(I_vals[i])
        curr_R_0 = 150.0 #FIXME
        if gamma == 0
            curr_R_0 == 150.0 #INF
        else
            curr_R_0 = beta/gamma
        end
        push!(beta_vals, beta)
        push!(gamma_vals, gamma)
        push!(R_0, curr_R_0)
    end
    #FIXME: last element is null
    push!(beta_vals, 0.0)
    push!(gamma_vals, 0.0)
    push!(R_0, 0.0)

    #add arrays to curr_df
    insertcols!(curr_df, 4, :cumulative_recovery=>cumulative_recovery)
    #curr_df.cumulative_death = cumulative_death
    curr_df.new_cases = new_cases
    curr_df.new_recovery = new_recovery
    curr_df.new_death = new_death
    curr_df.Infected = I_vals
    curr_df.Removed = R_vals
    curr_df.beta_t = beta_vals
    curr_df.gamma_t = gamma_vals
    curr_df.R_0 = R_0
    #@show curr_df


    #filepath = "/Users/evaxueyaoguo/Desktop/Florida/" * curr_county * ".csv"
    filepath = "./" * state_name * "/" * curr_county * ".csv"
    CSV.write(filepath, curr_df)

    #plot1
    p = plot(curr_df.date, beta_vals,
        title = curr_county * ", " * state_name * " β(t) and γ(t) from NYT data",
        label = "β(t)")
    plot!(curr_df.date, gamma_vals,
        label = "γ(t)")
    xlabel!("dates")
    savefig("./"*state_name*"/figs/1_" * curr_county * ".png")

    #plot2
    p = plot(curr_df.date, R_0,
        title = curr_county * ", " * state_name * " R0(t) from NYT data",
        label = "R0(t)")
    hline!([1.0])
    xlabel!("dates")

    savefig("./"*state_name*"/figs/2_" * curr_county * ".png")

end
