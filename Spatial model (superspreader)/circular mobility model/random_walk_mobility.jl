using Random, Distributions
using Plots
using CoordinateTransformations, StaticArrays

#Parameters
N = 1000 #150-900
L = 100.0 #m
r0 = 2.0
w0 = 1.0
λ = 0.2
#γ = 0.00046 #rate of recovery and death
tf = 60 #min #time unit unspecified (can be day or hour or even second) #40
model = "hub" #"hub" or "strong infectiousness"
fig = plot()
rw_interval = 1.0 #size of random walk steps

#--------------------------------------
function generate_random_event(p)
    @assert p >= 0 && p<=1
    u = rand(Float64)
    if u >= 0 && u < p #u:[0,p)-->event 1
        return 1
    else #u: [p,1)-->event 2
        return 0
    end
end

function dist(p1, p2)
    x1 = p1[1]
    x2 = p2[1]
    y1 = p1[2]
    y2 = p2[2]
    return (((x1-x2)^2) + ((y1-y2)^2))^0.5
end

function generate_infection_prob(r, α)
    #strong infectiousness model, w(r)=0 when r>r0
    if model == "strong infectiousness"
        if r <= r0
            w = w0 * ((1-(r/r0))^α)
        elseif r>r0
            w = 0
        end
    elseif model == "hub"
        if α == 2 #normal
            rn = r0
        elseif α == 0 #ss
            rn = (6.0^0.5)*r0
            α = 2
        end
        if r <= rn
            w = w0 * ((1-(r/rn))^α)
        elseif r > rn
            w = 0
        end
    end
    return w
end

function parse_pop(pop)
    S = Set(); I = Set(); R = Set()
    for i in 1:N
        if pop[i][4] == "S"
            push!(S, pop[i])
        elseif pop[i][4] == "I"
            push!(I, pop[i])
        elseif pop[i][4] == "R"
            push!(R, pop[i])
        end
    end
    return S, I, R
end

function initialize_pop(L, λ)
    @assert N>0
    pop = Any[] #might be faster to use a hash map
    #pss = generate_random_event(λ)
    pss = generate_random_event(1.0) # assume that the first individual is a superspreader
    if pss == 1 #ss
        p0 = [L/2, 0, "ss", "I", 1] #[x,y,"ss"/"n", "S"/"I"/"R", idx]
    elseif pss == 0 #normal
        p0 = [L/2, 0, "n", "I", 1]
    end
    push!(pop, p0)

    #generate the rest of the pop and their initial locations
    for i = 1:(N-1)
        x = rand(Float64)*L
        y = rand(Float64)*L
        pss = generate_random_event(λ)
        if pss ==1 #ss
            p = [x, y, "ss", "S", i+1] #Susceptible
        else #normal
            p = [x, y, "n", "S", i+1] #Susceptible
        end
        push!(pop, p)
    end
    return pop
end

function visualize(fig)
    #plot the location of each SIR individual
    s_x = Float64[]
    s_y = Float64[]
    i_x = Float64[]
    i_y = Float64[]
    r_x = Float64[]
    r_y = Float64[]
    ss_x = Float64[]
    ss_y = Float64[]
    for p in enumerate(S)
        curr = p[2]
        if length(curr)>8; x = curr[9]; y = curr[10]
        else x = curr[1]; y = curr[2]; end
        push!(s_x, x); push!(s_y, y)
    end
    for p in enumerate(I)
        curr = p[2]
        if length(curr)>8; x = curr[9]; y = curr[10]
        else x = curr[1]; y = curr[2]; end
        push!(i_x, x); push!(i_y, y)
    end
    for p in enumerate(R)
        curr = p[2]
        if length(curr)>8; x = curr[9]; y = curr[10]
        else x = curr[1]; y = curr[2]; end
        push!(r_x, x); push!(r_y, y)
    end

    pop = union(S, I, R)
    for p in enumerate(pop)
        curr = p[2]
        if curr[3]=="ss"
            if length(curr)>8; x = curr[9]; y = curr[10]
            else x = curr[1]; y = curr[2]; end
            push!(ss_x, x); push!(ss_y, y)
        end
    end

    plot!(fig, s_x, s_y, seriestype = :scatter, aspect_ratio = 1,
        #legend = false, xlims = (0.0, L), ylims = (0.0, L), label = "S")
        #legend = false, label = "S")
        legend = false, xlims = (0.0, L), ylims = (0.0, L), label = "S")
    plot!(fig, i_x, i_y, seriestype = :scatter, label = "I", c = :red)
    plot!(fig, r_x, r_y, seriestype = :scatter, label = "R", c = :yellow)
    plot!(fig, ss_x, ss_y, seriestype = :scatter, shape = :x, color = :black, label = "ss")
    title!("N="* string(N)* ", L="*string(L)*", λ="*string(λ)*
    ",\nRWstep="*string(rw_interval)*", tf="*string(tf))
    display(fig)
end

function infect() #modified based on ss_strong_infectiousness.jl
    for p in enumerate(I)
        for q in enumerate(S)
            r = dist(p[2], q[2])
            if p[2][3]=="ss"; α = 0 #ss
            else; α = 2 #normal
            end
            w = generate_infection_prob(r, α)
            event = generate_random_event(w) #1 = infection, 0 = nothing
            if event == 1 #infection
                curr = q[2]
                delete!(S, curr)
                push!(I, curr)
            end #else: nothing
        end
    end

    #=
    #recovery
    for p in enumerate(I)
        event = generate_random_event(γ)
        if event == 1 #recover/death
            curr = p[2]
            delete!(I, curr)
            push!(R, curr)
        end
    end
    =#
    #can add vaccination here

    println(length(S),", ", length(I),", ",length(R))
end

function move() #move 1 entire circle for now
    #for t in 1:time_steps #1 period/round of movement e.g. 1 day, 1 month etc.
        infect() #***
        fig = plot() #comment this line to see all locations through time
        visualize(fig)
        global S = rw_move_set(S)
        global I = rw_move_set(I)
        global R = rw_move_set(R)
        push!(S_vals, length(S)); push!(I_vals, length(I)); push!(R_vals, length(R))
    #end
end

function rw_move_set(set)
    set_new = deepcopy(set)
    for i in enumerate(set)
        p = i[2]
        x0 = p[1]; y0 = p[2] #consider boundary cases

        #random walk
        #x, y, +-1
        x_event = generate_random_event(0.5) #0 or 1 -->-1 or +1
        y_event = generate_random_event(0.5)

        x = 0.0
        y = 0.0
        if x_event == 0
            x = x0 - rw_interval
        elseif x_event==1
            x = x0 + rw_interval
        end

        if y_event == 0
            y = y0 - rw_interval
        elseif y_event==1
            y = y0 + rw_interval
        end

        #update x, y in set
        delete!(set_new, p)
        #println(length(set))
        #consider boundary cases
        if x < 0.0; x = 0.0
        elseif x > L; x = L; end
        if y < 0.0; y = 0.0
        elseif y > L; y = L; end
        #println("(", x, ", ", y, ")")
        @assert x>=0 && x<=L
        @assert y>=0 && y<=L
        p[1] = x; p[2] = y
        push!(set_new, p)
    end
    return set_new
end

#testing
pop = initialize_pop(L, λ)
S, I, R = parse_pop(pop)
S_vals = Float64[]; I_vals = Float64[]; R_vals = Float64[]

anim = @animate for i in 1:tf
    move() #e.g. 40 days
    #end of day: output
end
#gif(anim, "./figs/rw_2.gif", fps = 15)
