using Random, Distributions
using Plots
using CoordinateTransformations, StaticArrays
#FIXME: better plots

#Parameters
N = 10000 #150-900
loc_l = 100.0 #m
r0 = 2.0 #???
w0 = 1.0
λ = 0.3 #ratio of superspreaders
γ = 0.0 #rate of recovery and death
μse = 0.1 #FIXME: μse = β I(t)/N
ζ = .3 #E-->I
ρ = .03 #E-->L # probability of, once infected, going to L compartment
η = .8 #L-->ICU
κ = .3 #ICU-->R
tf = 5 #time unit unspecified (can be day or hour or even second) #40
model = "hub" #"hub" or "strong infectiousness"
time_steps = 12.0 #determines periodical movement
μ = 10.0
σ = 10.0
c_from_p = CartesianFromPolar()
p_from_c = PolarFromCartesian()
const global delta = (2 * π)/time_steps
hospital_x = loc_l * 1.5
hospital_y = loc_l * 0.5

#----------------------------
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

function circleShape(a, b, r)
    #https://discourse.julialang.org/t/plot-a-circle-with-a-given-radius-with-plots-jl/23295/15
    gr()
    θ = LinRange(0, 2*π, 500)
    a .+ r*sin.(θ), b .+ r*cos.(θ)
end

function parse_pop(pop)
    S = Set(); E = Set(); L = Set(); ICU = Set(); I = Set(); R = Set()
    for i in 1:N
        if pop[i][4] == "S"
            push!(S, pop[i])
        elseif pop[i][4] == "E"
            push!(E, pop[i])
        elseif pop[i][4] == "L"
            push!(L, pop[i])
        elseif pop[i][4] == "ICU"
            push!(ICU, pop[i])
        elseif pop[i][4] == "I"
            push!(I, pop[i])
        elseif pop[i][4] == "R"
            push!(R, pop[i])
        end
    end
    return S, E, L, ICU, I, R
end
#-------------------------------------------
function initialize_pop(loc_l, λ)
    @assert N>0
    pop = Any[] #might be faster to use a hash map
    pss = generate_random_event(λ)
    if pss == 1 #ss
        p0 = [loc_l/2, 0, "ss", "I"]
    elseif pss == 0 #normal
        p0 = [loc_l/2, 0, "n", "I"]
    end
    push!(pop, p0)

    #generate the rest of the pop and their initial locations
    for i = 1:(N-1)
        x = rand(Float64)*loc_l
        y = rand(Float64)*loc_l
        pss = generate_random_event(λ)
        if pss ==1 #ss
            p = [x, y, "ss", "S"] #Susceptible
        else #normal
            p = [x, y, "n", "S"] #Susceptible
        end
        push!(pop, p)
    end
    return pop
end

function visualize()
    p1 = plot()
    p2 = plot()
    #plot the location of each SIR individual
    s_x = Float64[]
    s_y = Float64[]
    i_x = Float64[]
    i_y = Float64[]
    r_x = Float64[]
    r_y = Float64[]
    e_x = Float64[]
    e_y = Float64[]
    l_x = Float64[]
    l_y = Float64[]
    icu_x = Float64[] #FIXME: better visualization for ICU
    icu_y = Float64[]
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

    for p in enumerate(E)
        curr = p[2]
        if length(curr)>8; x = curr[9]; y = curr[10]
        else x = curr[1]; y = curr[2]; end
        push!(e_x, x); push!(e_y, y)
    end

    for p in enumerate(L)
        curr = p[2]
        if length(curr)>8; x = curr[9]; y = curr[10]
        else x = curr[1]; y = curr[2]; end
        push!(l_x, x); push!(l_y, y)
    end

    for p in enumerate(ICU)
        curr = p[2]
        if length(curr)>8; x = curr[9]; y = curr[10]
        else x = curr[1]; y = curr[2]; end
        push!(icu_x, x); push!(icu_y, y)
    end

    pop1 = union(S, I, R); pop2 = union(E, L, ICU)
    pop = union(pop1, pop2)
    @assert length(pop)==N

    for p in enumerate(pop)
        curr = p[2]
        if curr[3]=="ss"
            if length(curr)>8; x = curr[9]; y = curr[10]
            else x = curr[1]; y = curr[2]; end
            push!(ss_x, x); push!(ss_y, y)
        end
        r = curr[6]
        a = curr[7]
        b = curr[8]
        ##plot!(p1, circleShape(a,b,r), seriestype = [:shape], lw = 0.5, c = :blue, linecolor = :black, legend = false, fillalpha = 0.2, aspect_ratio = 1)
        #plot!(p1, circleShape(a,b,r), seriestype = [:shape], lw = 0.5,c = :blue, linecolor = :black, fillalpha = 0.2, aspect_ratio = 1)
    end

    plot!(p1, s_x, s_y, seriestype = :scatter, aspect_ratio = 1,
        #legend = false, xlims = (0.0, loc_l), ylims = (0.0, loc_l), label = "S")
        #legend = false, label = "S")
        legend = false, xlims = (0.0, hospital_x), ylims = (0.0, loc_l), label = "S", c=:blue)
    plot!(p1, e_x, e_y, seriestype = :scatter, label = "E", c = :yellow)
    plot!(p1, l_x, l_y, seriestype = :scatter, label = "L", c = :orange)
    plot!(p1, icu_x, icu_y, seriestype = :scatter, label = "ICU", c = :purple)
    plot!(p1, i_x, i_y, seriestype = :scatter, label = "I", c = :red)
    plot!(p1, r_x, r_y, seriestype = :scatter, label = "R", c = :green)
    plot!(p1, ss_x, ss_y, seriestype = :scatter, shape = :x, color = :black, label = "ss")
    icu_annotate = "ICU" * string(length(ICU))
    s_annotate = "S " * string(length(S))
    e_annotate = "E " * string(length(E))
    l_annotate = "L " * string(length(L))
    i_annotate = "I " * string(length(I))
    r_annotate = "R " * string(length(R))
    n_annotate = "N " * string(length(pop))

    annotate!(0.9*hospital_x, 1.0*loc_l, n_annotate)
    annotate!(0.9*hospital_x, 0.9*loc_l, s_annotate)
    annotate!(0.9*hospital_x, 0.8*loc_l, e_annotate)
    annotate!(0.9*hospital_x, 0.7*loc_l, i_annotate)
    annotate!(0.9*hospital_x, 0.6*loc_l, l_annotate)
    annotate!(0.9*hospital_x, 0.5*loc_l, icu_annotate)
    annotate!(0.9*hospital_x, 0.4*loc_l, r_annotate)


    plot!(p2, S_vals, xlims = (0.0, tf*time_steps), ylims = (0.0, N)
        ,label = "S", legend=true)
    plot!(p2, I_vals, label = "I")
    plot!(p2, R_vals, label = "R")

    fig = plot(p1, p2, layout = (2,1))
    display(fig)

end

function generate_random_radii(μ, σ, N)
    #polling radius from a trcucated normal distribution
    #TODO: can use different distributions, or dist estimated from data
    truncated_normal = truncated(Normal(μ, σ), 0.0, Inf) #includive
    radii = rand(truncated_normal, N)
    for i in 1:N
        push!(pop[i], radii[i])
    end
end

function generate_random_vector(N)
    vectors = rand(Uniform(0, 2π), N)#[0, 2π]
    for i in 1:N
        push!(pop[i], vectors[i])
    end
end

function genetate_circular_path() #(radius, vector)-->(a,b)
    for i in 1:N
        p = pop[i]; x0 = p[1]; y0 = p[2]; θ = p[5]; r = p[6]
        dx = r * sin(θ)
        dy = r * cos(θ)
        #(a, b) is the center of the circular path
        #equation of the circular path: (x-a)^2 + (y-b)^2 = r^2
        a = x0 + dx
        b = y0 + dy
        push!(pop[i], a); push!(pop[i], b)
    end
end

function infect() #modified based on ss_strong_infectiousness.jl
    Infectious = union(I, L)
    #print(length(Infectious))
    for p in enumerate(Infectious)
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
                push!(E, curr)
            end #else: nothing
        end
    end

    for p in enumerate(E) #E-->L
        event = generate_random_event(ρ)
        if event == 1
            curr = p[2]
            delete!(E, curr)
            push!(L, curr)
        end
    end

    for p in enumerate(E) #E-->I
        event = generate_random_event(ζ)
        if event == 1
            curr = p[2]
            delete!(E, curr)
            push!(I, curr)
        end
    end

    for p in enumerate(L) #L-->ICU
        event = generate_random_event(η)
        if event == 1
            curr = p[2]
            delete!(L, curr)
            push!(ICU, curr)
        end
    end

    for p in enumerate(ICU) #ICU-->R
        event = generate_random_event(κ)
        if event == 1
            curr = p[2]
            delete!(ICU, curr)
            push!(R, curr)
        end
    end

    #recovery
    for p in enumerate(I) #I-->R
        event = generate_random_event(γ)
        if event == 1 #recover/death
            curr = p[2]
            delete!(I, curr)
            push!(R, curr)
        end
    end
    N = length(S) + length(E) + length(L) + length(ICU) + length(I) + length(R)
    println("N=", N,", S=", length(S),", E=", length(E), ", L=", length(L),
    ", ICU=", length(ICU), ", I=", length(I),", R=",length(R))
end


function move() #move 1 entire circle for now
    for t in 1:time_steps
        infect() #***
        visualize()

        global S = move_set(S)
        global E = move_set(E)
        global L = move_set(L)
        global ICU = move_set_ICU(ICU)
        global I = move_set(I)
        global R = move_set(R)

        push!(S_vals, length(S)); push!(I_vals, length(I)); push!(R_vals, length(R))
        push!(E_vals, length(E)); push!(L_vals, length(L)); push!(ICU_vals, length(ICU))
    end
end

function move_set(set)
    set_new = deepcopy(set)
    for i in enumerate(set)
        p = i[2]
        a = p[7]; b = p[8]; r = p[6]; x0 = p[1]; y0 = p[2] #consider the underlying location for boundary cases
        x_prime = x0 - a; y_prime = y0 - b

        xy = SVector(x_prime, y_prime)
        rθ = p_from_c(xy)
        @assert rθ.r ≈ r
        θ0 = rθ.θ

        rθ_new = Polar(r, θ0 + delta)
        xy_new = c_from_p(rθ_new)
        x = xy_new[1] + a; y = xy_new[2] + b

        #update x, y in set
        delete!(set_new, p)
        #println(length(set))
        p[1] = x; p[2] = y #underlying location
        #consider boundary cases
        if length(p)>8; pop!(p); pop!(p); end
        if x < 0.0; x = 0.0
        elseif x > loc_l; x = loc_l; end
        if y < 0; y = 0.0
        elseif y > loc_l; y = loc_l; end

        push!(p, x); push!(p, y) #visualized location (for boundary cases)
        push!(set_new, p)
    end
    return set_new
end

function move_set_ICU(ICU)
    set_new = deepcopy(ICU)
    for i in enumerate(ICU)
        p = i[2]
        a = p[7]; b = p[8]; r = p[6]; x0 = p[1]; y0 = p[2] #consider the underlying location for boundary cases
        delete!(set_new, p)
        if length(p)>8
            p[9]= hospital_x; p[10]=hospital_y
        else
            push!(p, hospital_x); push!(p, hospital_y)
        end
        push!(set_new, p)
    end
    return set_new
end

#testing
pop = initialize_pop(loc_l, λ)
generate_random_vector(N)
generate_random_radii(μ, σ, N)
genetate_circular_path()
S, E, L, ICU, I, R = parse_pop(pop)
S_vals = Float64[]; I_vals = Float64[]; R_vals = Float64[]
E_vals = Float64[]; L_vals = Float64[]; ICU_vals = Float64[]


@gif for i in 1:tf
    move() #e.g. 40 days
    #end of day: output
end

#can add herd immunity: vaccination
