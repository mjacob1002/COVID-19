using Random
#randomize the states of each state-change
function p_simulation(μ)
    u = rand()
    if 0<=u<=μ
        return true
    else
        return false
    end
end




# time is in days, updates the S, I, R every discrete day
function update_SIR(S0, I0, R0, V0, mu_SI, mu_IR, mu_RS, mu_SV, tm)
    S = zeros(Float64, tm)
    I = zeros(Float64, tm)
    R = zeros(Float64, tm)
    V = zeros(Float64, tm)
    t = zeros(Int64, tm)
    S[1] , I[1], R[1], V[1], t[1] = S0, I0, R0, 0,  0

    for i in 2:tm
        t[i] = i-1
        S_I = 0
        I_R = 0
        R_S = 0
        S_V = 0
        # state changes for S
        for j in 1:S[i-1]
            if(p_simulation(mu_SI) == true)
                S_I +=1
            elseif(p_simulation(mu_SV) == true)
                S_V +=1
            end
        end
        for j in 1:I[i-1]
            if(p_simulation(mu_IR) == true)
                I_R +=1
            end
        end
        for j in 1:R[i-1]
            if(p_simulation(mu_RS) == true)
                R_S +=1
            end
        end
        S[i] = S[i-1] + R_S - S_I - S_V
        I[i] = I[i-1] + S_I - I_R
        R[i] = R[i-1] + I_R - R_S
        V[i] = V[i-1] + S_V
    end
    return S, I, R, V, t
end

S0 = 99
I0 = 1
R0 = 0
V0 = 0
mu_SI = .1
mu_IR = .3
mu_RS = .25
mu_SV = .5
S, I, R, V, t = update_SIR(S0, I0, R0, V0, mu_SI, mu_IR, mu_RS, mu_SV, days)
ϵ = 1*10^(-12)

using Plots

 p = plot(t, S, label = "Susceptibles")
 plot!(t,I, label = "Infected")
 plot!(t,R, label = "Removed")
 plot!(t, V, label = "Vaccinated")
 xlabel!("Days")
 ylabel!("Number of People")
