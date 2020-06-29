# creates a 4 by 4 transition matrix for the standard DTMC SIR model
function trans_matrix(mu_SI, mu_IR, mu_RS, mu_SV)
    # S, I, R , V by S, I, R, V
    # p - transition matrix
    p = zeros(Float64, 4, 4)
    p[1,1] = (1-mu_SI-mu_SV)
    p[1,2] = (mu_SI)
    p[1,4] = (mu_SV)
    p[2,2] = (1-mu_IR)
    p[2,3] = (mu_IR)
    p[3,1] = mu_RS
    p[3,3] = 1- mu_RS
    p[4,4] = 1
    return p
end

# time is in days, updates the S, I, R every discrete day
function update_SIR(S0, I0, R0, V0, prob, tm)
    S = zeros(Float64, tm)
    I = zeros(Float64, tm)
    R = zeros(Float64, tm)
    V = zeros(Float64, tm)
    t = zeros(Int64, tm)
    S[1] , I[1], R[1], V[1], t[1] = S0, I0, R0, 0,  1

    for i in 2:tm
        t[i] = i
        S[i] = S[i-1] * prob[1,1] + R[i-1]*prob[3,1]
        I[i] = S[i-1]*prob[1,2] + I[i-1]*prob[2,2]
        R[i] = I[i-1] * prob[2,3] + R[i-1]*prob[3,3]
        V[i] = V[i-1] + S[i-1] * prob[1,4]
    end
    return S, I, R, V, t
end

S0 = 99
I0 = 1
R0 = 0
V0 = 0
mu_SI = .10
mu_IR = .2
mu_RS = .05
mu_SV = .3
days = 30
m = trans_matrix(mu_SI, mu_IR, mu_RS, mu_SV)
S, I, R, V, t = update_SIR(S0, I0, R0, V0, m, days)
ϵ = 1*10^(-12)

using Plots

 p = plot(t, S, label = "Susceptibles")
 plot!(t,I, label = "Infected")
 plot!(t,R, label = "Removed")
 plot!(t, V, label = "Vaccinated")
 xlabel!("Days")
 ylabel!("Number of People")
