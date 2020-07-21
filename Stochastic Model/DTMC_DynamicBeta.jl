using Plots
using Random
# other variables
global S_0 = 1000000 * 8.4
global I_0 =1
global R_0 = 0
global B0 = 5.5
#starting hospital capacity at time = 0
global C0 = .003 *S_0
global N = S_0 + I_0 + R_0


function p_simulation(μ,t)
     x = rand()
    if(0<=x<=μ)
        return true
    else
        return false
    end
end
# will act as multiplier to probability to vary the chance based on
#current conditions
function my_alpha(t, I)
    capacity = hosp_limit(C0,t)
    local m = I/capacity
    α = 1
    if(.35< m < .5)
        α = .65
    elseif .5 <= m < .8
        α = .35
    elseif .8 <= m
        α = .05
    end
    return α
end

function update_SIR(time, mu_IR)
    num_int = time + 1
    # initialize the arraysS
    S = zeros(Float64, num_int)
    I = zeros(Float64, num_int)
    R = zeros(Float64, num_int)
    t = zeros(Int64, num_int)
    S[1] = S_0
    I[1] = I_0
    R[1] = R_0
    t[1] = 0
    for i in 2 : num_int
        t[i] = i-1
        mu_SI = my_muSI(t[i-1], I[i-1])
        SI = 0
        IR = 0
        print("mu_SI: ", mu_SI)
        print(" alpha: ", my_alpha(i-1, I[i-1]))
        println( " ratio: ", I[i-1]/C0)
        for j in 2:S[i-1]
            if(p_simulation(mu_SI, i-1))
                SI +=1
            end
        end
        for j in 2:I[i-1]
            if(p_simulation(mu_IR, i-1))
                IR+=1
            end
        end
        S[i] = S[i-1] - SI
        I[i] = I[i-1] + SI - IR
        R[i] = R[i-1] + IR
    end
    return t, S, I , R
end

#accounts for the constantly increasing hospital limit as function of time
#currently a linear function, will experiment with other functions
function hosp_limit(C_0, t)
    c = C0 + .00075*t
    return c
end

function my_muSI(t, I)
    α = my_alpha(t, I)
    return α*B0* I/N
end
#basic SIR model(as a DTMC)
function main()
    mu_IR = 3/15
    duration = 365
    t, S, I, R = update_SIR(duration, mu_IR)
    p = plot(t, I, label = "Infected")
    title!("Dynamic DTMCv3")
end
println("New Run: ")
main()
