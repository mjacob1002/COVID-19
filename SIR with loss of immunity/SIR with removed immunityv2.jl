using StatsPlots
using Plots
# returns derivatives as a vector containing dS/dt,dI/dt,dR/dt
#1/v is the typical time it takes to lose immunity
function f(β, γ, u, t, v)
    S, I, R = u
    N = S + I + R
    return [-β*S*I/N + v*R,β*S*I-γ*I,γ*I-v*R]
end
# eventually will make beta a function of other variables for "the dance"
function beta()
    return
end

function R_naught(β, τ)
    return β*τ
end

function ode_FE(u_0, dt, T, Β, γ, v)
    # the number of time intervals(T/dt)
    N_t = convert(Int64, round(T/dt)) # integer amount of intervals
    nrows = N_t+1 # because the indexing starts at 1
    ncolumns = length(u_0)
    # creates a 2D array to store all of the S,I,R values
    u = zeros(Float64, nrows, ncolumns)
    # time frame to be plotted
    t = range(0, N_t*dt, length = nrows)
    #initialize the first row with the starting values
    u[1,:] = u_0
    #Stepping forward using Euler's Method
    for i in 1:N_t
        u[i+1,:] = u[i,:] + dt*f(Β, γ, u[i,:], t[i], v)
    end
    #returns the series with completed S,I,R, as well as t
    return u,t
end

function demo()
#Initial S, I, R
    S0 = 50
    I0 = 1
    R0 = 0
    U_0 = [S0, I0, R0]
    #parameters
    Β = 10.0/(40*8*4) # time frame of 1 day; P(infected and meet in pair)
    γ = 3.0/15
    dt = 1/24 # 1 hour
    D = 1200 #Numbers of days being simulated
    v = 1/100
    u,t = ode_FE(U_0, dt, D, Β, γ, v)
    #moving the S,I,R from the big matrix to vectors S,I,R
    S = u[: , 1]
    I = u[:, 2]
    R = u[:,3]
    # plotting the graph
    p = plot(t,S, label = "Susceptibles")
    plot!(t,I, label = "Infected")
    plot!(t,R, label = "Removed")
    xlabel!("Daysv2")
    ylabel!("Number of People")
    title!("Generalized Spreading of Disease")
    #=Consistency check
    N = S0 + I0 + R0
    ϵ = 1E-12
    for i in 1:length(S)
        SIR_sum = S[i] + I[i] + R[i]
        if abs(SIR_sum - N) > ϵ
            println("Consistency Check Failed: S+I+R = %g != %g", (SIR_sum, N))
        end
    end
    =#
end

demo()
