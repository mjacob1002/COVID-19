using DataFrames, CSV
using Statistics
using Plots

lamda_0 = [2.2, 3.44, 4.61, 6.24, 8.5, 10.9, 12.14, 13.0, 13.96, 20.93]
lamda_2 = [2.02, 2.9, 3.35, 4.18, 5.16, 6.84, 9.61, 11.82, 13.04, 20.39]
lamda_4 = [2.0, 2.68, 3.05, 3.65, 4.33, 5.43, 7.25, 10.46, 12.68, 19.95]
lamda_6 = [2.0, 2.26, 3.0, 3.16, 3.97, 4.7, 5.94, 8.67, 12.12, 19.32]
lamda_8 = [2.0, 2.07, 2.97, 3.04, 3.74, 4.23, 5.29, 7.44, 11.47, 18.95]
lamda_10 = [2.0, 2.06, 2.85, 3.0, 3.35, 4.07, 4.97, 6.7, 10.96, 19.01]

x = [10.0:10.0:100.0]

fig = plot(x, lamda_0, label = "λ=0.0", xticks=10.0:10.0:100.0, legend = :topleft,
    )
plot!(fig, x, lamda_2, label = "λ=0.2")
plot!(fig, x, lamda_4, label = "λ=0.4")
plot!(fig, x, lamda_6, label = "λ=0.6")
plot!(fig, x, lamda_8, label = "λ=0.8")
plot!(fig, x, lamda_10, label = "λ=1.0")
xlabel!(fig, "%ppl infected (cumulative)")
ylabel!(fig, "time")
title!("Time vs %ppl infected \nfor different superspreader ratios")
savefig("/Users/xueyaoguo/Desktop/Monte_Carlo/time_percent_3.png")
