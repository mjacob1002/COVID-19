using StatsPlots
using Plots

# returns the derivative as a vector
function f(u,t, β, γ)
    S, I, R = u
    return [-β*S*I, β*S*I-γ*I, γ*I]
end

function demo()

    println("Demo ")
    β = 10.0 / (5*40*24)
    γ = 3.0 / (15*24)
    dt = .1 #  6 mins(.1 of an hour)
    D = 5 # Simulate for D days
    N_t = convert(Int64,(24D/dt)) #number of points for model
    T = dt*N_t
    S0 = 50
    I0 = 1
    R0 = 0
    U_0 = [S0, I0, R0]

    println("N_t: ", N_t)
    println("dt: ", dt)
    println("T: ", T)

    u,t = ode_FE(U_0, dt, T, β, γ)

    S = u[:, 1]
    I = u[:, 2]
    R = u[:, 3]

    p = plot(t,S, label = "Susceptible")
    plot!(t,I, label = "Infected")
    plot!(t,R, label = "Removed")
    xlabel!("Hours")
    ylabel!("Number of People")
    title!("Spreading of Disease")

    #Consistency check
    N = S0 + I0 + R0
    ϵ = 1E-6
    for i in 1:length(S)
        SIR_sum = S[i] + I[i] + R[i]
        if abs(SIR_sum - N) > ϵ
            println("Consistency Check Failed: S+I+R = %g != %g", (SIR_sum, N))
        end
    end
end

# Euler's Method
function ode_FE(u_0, dt, T, β, γ)
    # number of time intervals [ = T/dt]
    N_t = convert(Int64, round(convert(Float64, T/dt)))

    #create 2D array to hold all values
    nrows = N_t+1
    ncols = length(u_0)
    u = zeros(Float64, nrows, ncols)

    t= range(1, N_t*dt,length=nrows)

    u[1,:] = u_0
    for i in 1:N_t
        println("i: ", i)
        r_0 = f(u[i,:], t[i], β, γ)
        r_0 = dt*r_0
        u[i+1,:] = u[i,:] + r_0
    end

    return u,t
end


demo()
