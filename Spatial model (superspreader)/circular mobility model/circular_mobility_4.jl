using Random, Distributions
using Plots

#Parameters
N = 10
L = 100.0 #m
r0 = 2.0
w0 = 1.0
λ = 0.3
γ = 0.1 #rate of recovery and death
tf = 40 #time unit unspecified (can be day or hour or even second)
model = "hub" #"hub" or "strong infectiousness"
time_steps = 12


function initialize_pop(L, λ)
    #FIXME: ASSERT N>0
    pop = Any[] #might be faster to use a hash map
    #I0 = 1
    pss = generate_random_event(λ)
    if pss == 1 #ss
        p0 = [L/2, 0, "ss", "I"]
    elseif pss == 0 #normal
        p0 = [L/2, 0, "n", "I"]
    end
    push!(pop, p0)
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
        push!(pop, p)
    end
    return pop
end

function visualize(pop, paths, fig) #tested
    #plot the location of each individual
    x_vals = Float64[]
    y_vals = Float64[]
    for i in enumerate(pop)
        #println(i[2])
        push!(x_vals, i[2][1])
        push!(y_vals, i[2][2])
    end

    for i in 1:N
        a = paths[i][1]
        b = paths[i][2]
        r = paths[i][3]
        #plot!(fig, circleShape(a,b,r), seriestype = [:shape], lw = 0.5, c = :blue, linecolor = :black, legend = false, fillalpha = 0.2, aspect_ratio = 1)
        plot!(fig, circleShape(a,b,r), seriestype = [:shape], lw = 0.5, linecolor = :black, legend = false, fillalpha = 0.2, aspect_ratio = 1)
    end

    plot!(fig, x_vals, y_vals, seriestype = :scatter, aspect_ratio = 1,
        legend = false, xlims = (0.0, L), ylims = (0.0, L))
    display(fig)

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

function generate_random_radii(μ, σ, N)
    #polling radius from a trcucated normal distribution
    #TODO: can use different distributions, or dist estimated from data
    truncated_normal = truncated(Normal(μ, σ), 0.0, Inf) #includive
    radii = rand(truncated_normal, N)
    return radii
end

function generate_random_vector(N)
    vectors = rand(Uniform(0, 2π), N)#[0, 2π]
    return vectors
end

function genetate_circular_path(pop, radii, vectors) #tested
    paths = Any[]
    for i in 1:N
        r = radii[i]
        θ = vectors[i]
        dx = r * sin(θ)
        dy = r * cos(θ)
        x0 = pop[i][1]
        y0 = pop[i][2]
        #(a, b) is the center of the circular path
        #equation of the circular path: (x-a)^2 + (y-b)^2 = r^2
        a = x0 + dx
        b = y0 + dy
        push!(paths, [a, b, r])
    end
    return paths
end

function circleShape(a, b, r)
    #https://discourse.julialang.org/t/plot-a-circle-with-a-given-radius-with-plots-jl/23295/15
    gr()
    θ = LinRange(0, 2*π, 500)
    a .+ r*sin.(θ), b .+ r*cos.(θ)
end

function move(paths) #move 1 entire circle for now
    #update locations in pop
    delta = (2 * π)/time_steps #dθ
    for t in 1:time_steps
        for i in 1:N
            #use polar coordinates with origin at (a,b)
            a = paths[i][1]
            b = paths[i][2]
            r = paths[i][3]
            x0 = pop[i][1]
            y0 = pop[i][2]

            x_prime = x0 - a
            y_prime = y0 - b
            θ0 = atan(y_prime/x_prime) #arctan:[-π, π]
            #4 quadrants
            if x_prime>=0
            elseif x_prime<0
                θ0 = θ0 + π
            end

            x = r * cos(θ0 + delta) + a
            y = r * sin(θ0 + delta) + b

            global pop[i][1] = x
            global pop[i][2] = y
        end

        #=visualization and boundary cases
        #FIXME:Note boundary cases: set x=L, y=L
        vis_pop
        =#
        vis_pop = copy(pop)
        for j in 1:N
            x = pop[j][1]; y = pop[j][2]
            if x < 0; vis_pop[j][1] = 0; end
            if x > L; vis_pop[j][1] = L; end
            if y < 0; vis_pop[j][2] = 0; end
            if y > L; vis_pop[j][2] = L; end
        end
        fig = plot()
        visualize(vis_pop, paths, fig)
    end
end

#add infection here

#test
pop = initialize_pop(L, λ)
radii = generate_random_radii(10.0, 1.0, N)
vectors = generate_random_vector(N)
paths = genetate_circular_path(pop, radii, vectors)
fig = plot()
visualize(pop, paths, fig)
move(paths)
