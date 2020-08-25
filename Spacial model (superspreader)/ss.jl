using Random
using Plots

#Parameters
N = 200 #150-900
r0 = 2.0
L = 10 * r0
w0 = 1.0
γ = 0.0 #recovery or death
λ = 0.3 #ss density

Infected = Set()
Susceptible = Set()
Removed = Set()

function dist(p1, p2)
    x1 = p1[1]
    x2 = p2[1]
    y1 = p1[2]
    y2 = p2[2]
    return (((x1-x2)^2) + ((y1-y2)^2))^0.5
end

function generate_random_event(p)
    #FIXME: assert u is [0,1]
    u = rand(Float64)
    if u >= 0 && u < p #u:[0,p)-->event 1
        return 1
    else #u: [p,1)-->event 2
        return 0
    end
end

function generate_infection_prob(r, α)
    #strong infectiousness model, w(r)=0 when r>r0
    if r<= r0
        w = w0 * ((1-(r/r0))^α)
    elseif r>r0
        w = 0
    end
    return w
end

function initialize_run(N, L)
    #determine if p0 is ss
    pss = generate_random_event(λ)
    if pss ==1 #ss
        p0 = (L/2, 0.0, "ss") #Infected
    else #normal
        p0 = (L/2, 0.0, "n") #Infected
    end
    push!(Infected, p0)

    for i = 1:(N-1)
        x = rand(Float64)*L
        y = rand(Float64)*L
        pss = generate_random_event(λ)
        if pss ==1 #ss
            p = (x, y, "ss") #Susceptible
        else #normal
            p = (x, y, "n") #Susceptible
        end
        push!(Susceptible, p)
    end
end

#Simulation
initialize_run(N, L) # a set of all ppl

arrows = Any[]
function simulate()
    for p in enumerate(Infected)
        for q in enumerate(Susceptible)
            r = dist(p[2], q[2]) #determine if p will infect q or not
            if p[2][3]=="ss"
                α = 0 #ss
            else
                α = 2 #normal
            end
            w = generate_infection_prob(r, α)
            event = generate_random_event(w) #1 = infection, 0 = nothing
            if event == 1 #infection
                curr = q[2]
                delete!(Susceptible, curr)
                push!(Infected, curr)
                #record arrow from p to q
                arr = ([p[2][1], q[2][1]],[p[2][2], q[2][2]])
                push!(arrows, arr)
            else #nothing
            end
        end
    end

    #recovery
    for p in enumerate(Infected)
        event = generate_random_event(γ)
        if event == 1 #recover/death
            curr = p[2]
            delete!(Infected, curr)
            push!(Removed, curr)
        end
    end

    #can add vaccination here
end

#simulation ends when
#1. percolation
#2. time_max = 40 ***
#3. infection = 0

function visualize()
    #Visualization
    s_x = Float64[]
    s_y = Float64[]
    i_x = Float64[]
    i_y = Float64[]
    r_x = Float64[]
    r_y = Float64[]
    for p in enumerate(Susceptible)
        curr = p[2]
        x = curr[1]
        y = curr[2]
        push!(s_x, x)
        push!(s_y, y)
    end
    for p in enumerate(Infected)
        curr = p[2]
        x = curr[1]
        y = curr[2]
        push!(i_x, x)
        push!(i_y, y)
    end
    for p in enumerate(Removed)
        curr = p[2]
        x = curr[1]
        y = curr[2]
        push!(r_x, x)
        push!(r_y, y)
    end

    p = plot(s_x, s_y, seriestype = :scatter, label = "S")
    plot!(p, i_x, i_y, seriestype = :scatter, label = "I")
    plot!(p, r_x, r_y, seriestype = :scatter, label = "R")

    for a in enumerate(arrows)
        Plots.display(plot!(a[2][1], a[2][2], label = false, linecolor = :red))
    end


    savefig("./figs/run_single.png")
end

for i in 1:40
    simulate()
end
visualize()
