using Random, Distributions
using Plots
using CoordinateTransformations, StaticArrays

function generate_random_event(p)
    #=
    function to generate a random event using a random number generator
    @param p: the probability of an event
    @return: 1 if the event happens, and 0 if the event does not happen
    =#
    @assert p >= 0 && p<=1
    u = rand(Float64)
    if u >= 0 && u < p #u:[0,p)-->event happens with probability p
        return 1
        else #u: [p,1)-->event does not happen with probability 1-p
        return 0
    end
end

function dist(p1, p2)
    #=
    function to get the distance between two people
    @param p1: the first person as an array of [x1, y1, ...]
    @param p1: the second person as an array of [x2, y2, ...]
    @return: the distance between p1 and p2
    =#
    x1 = p1[1] #get the coordinates for p1 and p2
    x2 = p2[1]
    y1 = p1[2]
    y2 = p2[2]
    return (((x1-x2)^2) + ((y1-y2)^2))^0.5 #distance 2d
end

function generate_infection_prob(r, α)
    #=
    Function to get the infection probability based on the two models in Fujie&Odagaki (2007)
    @param r: distance between the two people
    @param α: power param for the strong infectiousness model (if hub model, α=2.0)
    @return: infection probability calculated based on the infection probability in Fujie&Odagaki (2007)
    =#
    if model == "strong infectiousness" #strong infectiousness model, w(r)=0 when r>r0
        if r <= r0
            w = w0 * ((1-(r/r0))^α)
        elseif r>r0
            w = 0
        end
    elseif model == "hub" #hub model, α=2.0
        #for a detailed description of this implementation, see Fujie&Odagaki (2007) pp. 845
        if α == 2 #normal (α as a marker for normal spreader vs ss here, and α is changed to 2 at line 52 for the hub model)
            rn = r0
        elseif α == 0 #ss
            rn = (6.0^0.5)*r0
            α = 2 #set α to 2
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
    #=
    Function to parse an initialized population into SIR sets (or with addition sets like E, L, and ICU)
    @param pop: a population (an array of individual arrays) to be parsed
    @return: S, E, L, ICU, I, R sets containing their corresponding individuals
    =#
    S = Set(); E = Set(); L = Set(); ICU = Set(); I = Set(); R = Set() #initialize sets
    for i in 1:N #for every individual in the population
        if pop[i][4] == "S" #parse into sets following the set marker (pop[i][4])
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

function initialize_pop(loc_l, λ, N)
    #=
    Function to initialze a simulated population
    @param loc_l: size for the parameter space
    @param λ: superspreader ratio
    @param N: population size
    @return: an array of individuals representing the population
    =#
    @assert N>0
    pop = Any[] #might be faster to use a hash map
    #pss = generate_random_event(λ) #no assumptions about whether patient0 is a superspreader or not
    pss = generate_random_event(1.0) #assume patient0 is a superspreader
    if pss == 1 #ss
        p0 = [loc_l/2, 0, "ss", "I"] #an array representing an individual, can be considered a person object
    elseif pss == 0 #normal
        p0 = [loc_l/2, 0, "n", "I"]
    end
    push!(pop, p0) # add patient0 to population

    #generate the rest of the pop (susceptible)
    for i = 1:(N-1)
        x = rand(Float64)*loc_l #x and y coordinates are randomly generated within the parameter space
        y = rand(Float64)*loc_l
        pss = generate_random_event(λ)#determine if this person is ss or not
        if pss ==1 #ss
            p = [x, y, "ss", "S"] #mark as susceptible
        else #normal
            p = [x, y, "n", "S"] 
        end
        push!(pop, p) #add person generated to population
    end
    return pop
end

function generate_random_radii(pop, μ, σ, N)
    #=
    Function to generate random radii for the periodic pattern (circular movement in parameter space)
    The radii are pulled from a trcucated normal distribution.
    @TODO: can use different distributions, or dist estimated from data
    @param pop: the population that the random radii are generated for 
    @param μ: mean of the normal distribution
    @param σ: standard deviation of the normal distribution
    @param N: population size
    @return: the population with radii generated for the circular movement pattern
    =#
    truncated_normal = truncated(Normal(μ, σ), 0.0, Inf) #inclusive
    radii = rand(truncated_normal, N)
    for i in 1:N
        push!(pop[i], radii[i]) #add the radius to each person in the population
    end
    return pop
end

function generate_random_vector(pop, N)
    #=
    Function to generate random vectors for the periodic pattern (i.e. determining the starting location of the person on the circle)
    @param pop: the population that the vectors are generated for
    @param N: population size
    @return: the population with vectors/starting locations generated for the circular movement pattern
    =#
    vectors = rand(Uniform(0, 2π), N)#pulling uniformly from [0, 2π], each angle represents a starting location on the circle (aka a vector)
    for i in 1:N
        push!(pop[i], vectors[i])#add the vector to each person in the population
    end
    return pop
end

function genetate_circular_path(pop) #tested: (radius, vector)-->(a,b)
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
    return pop
end

function infect(S, E, L, ICU, I, R, N) #modified based on ss_strong_infectiousness.jl
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

    @assert length(S) + length(E) + length(L) + length(ICU) + length(I) + length(R) == N
    #=
    println("N=", N,", S=", length(S),", E=", length(E), ", L=", length(L),
    ", ICU=", length(ICU), ", I=", length(I),", R=",length(R))
    =#
    return S, E, L, ICU, I, R
end

function move(S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals, S, E, L, ICU, I, R, N) #move 1 entire circle for now #FIXME
    for t in 1:time_steps #1 period/round of movement e.g. 1 day, 1 month etc.
        S, E, L, ICU, I, R = infect(S, E, L, ICU, I, R, N) #***
        #global counter = counter + 1
        #visualize(counter)

        S = move_set(S)
        E = move_set(E)
        L = move_set(L)
        ICU = move_set_ICU(ICU)
        I = move_set(I)
        R = move_set(R)

        push!(S_vals, length(S)); push!(I_vals, length(I)); push!(R_vals, length(R))
        push!(E_vals, length(E)); push!(L_vals, length(L)); push!(ICU_vals, length(ICU))

    end

    return S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals, S, E, L, ICU, I, R
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