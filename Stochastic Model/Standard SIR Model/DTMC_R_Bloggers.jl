using Random
# adds the randomness required in markov chain
function p_sim(μ)
    x = rand()
    if 0<=x<=μ
        return true
    else
        return false
    end
end

# time is in days, updates the S, I, R every discrete day
function update_SIR(S0, I0, R0, mu_SI, mu_IR, tm)
    # create the matrix
    S = zeros(Float64, tm)
    I = zeros(Float64, tm)
    R = zeros(Float64, tm)
    t = zeros(Int64, tm)
    # initialize the values of the array
    S[1] , I[1], R[1], t[1] = S0, I0, R0, 0

    for i in 2:tm
        t[i] = i-1
        # keeps track on the change from S -> I, I -> R
        S_I = 0
        I_R = 0
        for j in 1:S[i-1]
            if(p_sim(mu_SI) == true)
                S_I +=1
            end
        end
        # adjusts the next S value for the next discrete day
        S[i] = S[i-1] - S_I
        for j in 1:I[i-1]
            if(p_sim(mu_IR) == true)
                I_R +=1
            end
        end
        # adjusts the I, R values for the next discrete day
        I[i] = I[i-1] + S_I - I_R
        R[i] = R[i-1] + I_R
    end
    return S, I, R, t
end

S0 = 99
I0 = 1
R0 = 0
mu_SI = .10
mu_IR = .2
days = 31
S, I, R, t = update_SIR(S0, I0, R0, mu_SI, mu_IR, days)


using Plots

 p = plot(t, S, label = "Susceptibles")
 plot!(t,I, label = "Infected")
 plot!(t,R, label = "Removed")
 xlabel!("Daysv2")
 ylabel!("Number of People")
