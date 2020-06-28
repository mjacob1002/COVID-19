# creates a 3 by 3 transition matrix for the standard DTMC SIR model
function trans_matrix(mu_SI, mu_IR)
    # S, I, R by S, I, R
    # p - transition matrix
    p = zeros(Float64, 3, 3)
    p[1,1] = 1-mu_SI
    p[1,2] = mu_SI
    p[3,3] = 1
    p[2,2] = 1 - mu_IR
    p[2,3] = mu_IR
    return p
end

# time is in days, updates the S, I, R every discrete day
function update_SIR(S0, I0, R0, prob, tm)
    S = zeros(Float64, tm)
    I = zeros(Float64, tm)
    R = zeros(Float64, tm)
    t = zeros(Int64, tm)
    S[1] , I[1], R[1], t[1] = S0, I0, R0, 1

    for i in 2:tm
        t[i] = i
        S[i] = S[i-1] * prob[1,1]
        I[i] = I[i-1] + prob[1,2]*S[i-1] - prob[2,3] * I[i-1]
        R[i] = R[i-1] + prob[2,3]*I[i-1]
    end
    return S, I, R, t
end

S0 = 99
I0 = 1
R0 = 0
mu_SI = .10
mu_IR = .2
days = 30
m = trans_matrix(mu_SI, mu_IR)
S, I, R, t = update_SIR(S0, I0, R0, m, days)


using Plots

 p = plot(t, S, label = "Susceptibles")
 plot!(t,I, label = "Infected")
 plot!(t,R, label = "Removed")
 xlabel!("Daysv2")
 ylabel!("Number of People")
