using Random, Distributions
using Plots
using CoordinateTransformations, StaticArrays

#Parameters
N = 1000 #150-900
L = 100.0 #m
r0 = 2.0 #???
w0 = 1.0
λ = 0.3
γ = 0.2 #rate of recovery and death
tf = 10 #time unit unspecified (can be day or hour or even second) #40
model = "hub" #"hub" or "strong infectiousness"
time_steps = 12.0 #determines periodical movement
μ = 3.0
σ = 1.0
fig = plot()
c_from_p = CartesianFromPolar()
p_from_c = PolarFromCartesian()
const global delta = (2 * π)/time_steps

#--------------unchanged functions----------
function generate_random_event(p)
    @assert p >= 0 && p<=1 #u is [0,1] #FIXME: assert
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
#-------------------------------------------
function initialize_pop(L, λ)
    @assert N>0 #FIXME: ASSERT N>0
    pop = Any[] #might be faster to use a hash map
    pss = generate_random_event(λ)
    if pss == 1 #ss
        p0 = [L/2, 0, "ss", "I"]
    elseif pss == 0 #normal
        p0 = [L/2, 0, "n", "I"]
    end
    push!(pop, p0) #; push!(Infected, p0)

    #generate the rest of the pop and their initial locations
    for i = 1:(N-1)
        x = rand(Float64)*L
        y = rand(Float64)*L
        pss = generate_random_event(λ)
        if pss ==1 #ss
            p = [x, y, "ss", "S"] #Susceptible
        else #normal
            p = [x, y, "n", "S"] #Susceptible
        end
        push!(pop, p) #; push!(Susceptible, p)
    end
    return pop
end

function visualize(fig) #tested
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
        r = curr[6]
        a = curr[7]
        b = curr[8]
        ##plot!(fig, circleShape(a,b,r), seriestype = [:shape], lw = 0.5, c = :blue, linecolor = :black, legend = false, fillalpha = 0.2, aspect_ratio = 1)
        #plot!(fig, circleShape(a,b,r), seriestype = [:shape], lw = 0.5,c = :blue, linecolor = :black, fillalpha = 0.2, aspect_ratio = 1)
    end

    plot!(fig, s_x, s_y, seriestype = :scatter, aspect_ratio = 1,
        #legend = false, xlims = (0.0, L), ylims = (0.0, L), label = "S")
        #legend = false, label = "S")
        legend = false, xlims = (0.0, L), ylims = (0.0, L), label = "S")
    plot!(fig, i_x, i_y, seriestype = :scatter, label = "I", c = :red)
    plot!(fig, r_x, r_y, seriestype = :scatter, label = "R", c = :yellow)
    plot!(fig, ss_x, ss_y, seriestype = :scatter, shape = :x, color = :black, label = "ss")
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

function genetate_circular_path() #tested: (radius, vector)-->(a,b)
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

    #recovery
    for p in enumerate(I)
        event = generate_random_event(γ)
        if event == 1 #recover/death
            curr = p[2]
            delete!(I, curr)
            push!(R, curr)
        end
    end
    #can add vaccination here

    println(length(S),", ", length(I),", ",length(R))
end


function move() #move 1 entire circle for now #FIXME
    for t in 1:time_steps #1 period/round of movement e.g. 1 day, 1 month etc.
        infect() #***
        fig = plot() #comment this line to see all locations through time
        visualize(fig)
        global S = move_set(S)
        global I = move_set(I)
        global R = move_set(R)
        push!(S_vals, length(S)); push!(I_vals, length(I)); push!(R_vals, length(R))
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
        if x < 0.0; x = 0.0; end
        if x > L; x = L; end
        if y < 0.0; y = 0.0; end
        if y > L; y = L; end
        push!(p, x); push!(p, y) #visualized location (for boundary cases)
        push!(set_new, p)
    end
    return set_new
end

#testing
pop = initialize_pop(L, λ)
generate_random_vector(N)
generate_random_radii(μ, σ, N)
genetate_circular_path()
S, I, R = parse_pop(pop)
S_vals = Float64[]; I_vals = Float64[]; R_vals = Float64[]

for i in 1:tf
    move() #e.g. 40 days
    #end of day: output
end

#can add herd immunity: vaccination
#exposed? and ICU
