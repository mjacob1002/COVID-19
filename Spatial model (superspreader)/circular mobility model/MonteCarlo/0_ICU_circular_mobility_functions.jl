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
    if u > 0 && u <= p #u:[0,p)-->event happens with probability p
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

function generate_infection_prob(r, α, r0)
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

function genetate_circular_path(pop)
    #=
    Function to get the center of the circular path: (radius, vector)-->(a,b)
    @param pop: a population with the radii and vectors generated
    @return: the population with the center of the circles stored in each individual array
    =#
    for i in 1:N
        p = pop[i]; x0 = p[1]; y0 = p[2]; θ = p[5]; r = p[6] #get parameters
        #θ is from generate_random_vector and is treated as the angle between the y axis and the line connecting (x0, y0) and (a,b)
        dx = r * cos(θ)
        dy = r * sin(θ)
        #(a, b) is the center of the circular path
        #equation of the circular path: (x-a)^2 + (y-b)^2 = r^2
        a = x0 - dx
        b = y0 - dy
        push!(pop[i], a); push!(pop[i], b)
    end
    return pop
end

function infect(S, E, L, ICU, I, R, N, r0 ,ρ, ζ, η, κ, γ)
    #=
    Function to update the population with infection, recovery, and other events in the compartmental model
    @params: sets to be updated
    @return: sets after the update
    =#
    Scopy = deepcopy(S)
    Ecopy = deepcopy(E)
    Lcopy = deepcopy(L)
    ICUcopy = deepcopy(ICU)
    Icopy = deepcopy(I)
    Rcopy = deepcopy(R)
    Infectious = union(I, L) #All infectious individuals = I + L
    #S-->E (an exposure event)
    for p in enumerate(Infectious) #for every infectious person
        for q in enumerate(S) #see if each susceptible person is exposed
            r = dist(p[2], q[2]) #get the distance between the two people
            if p[2][3]=="ss"; α = 0 #ss #FIXME: use a better marker
            else; α = 2 #normal
            end
            w = generate_infection_prob(r, α, r0) #get the infection probability (exposure probability in this case)
            event = generate_random_event(w) #1 = infection, 0 = nothing
            if event == 1 #infection
                curr = q[2]
                delete!(Scopy, curr) #the susceptible person is exposed and removed from S
                push!(Ecopy, curr) #and added to E
            end #else: nothing
        end
    end

    for p in enumerate(E) #E-->Infectious (infection that will go to the ICU)
        event = generate_random_event(ζ) #probablity = probability of leaving E
        if event == 1 # leaves E
            curr = p[2]
            delete!(Ecopy, curr)
            event2 = generate_random_event(ρ) # test whether goes to L or I
            if event2 == 1 # goes to L
                push!(Lcopy, curr)
            else # goes to I
                push!(Icopy, curr)

            end
        end
    end
    #=
    for p in enumerate(E) #E-->I (infection that won't go to the ICU)
        event = generate_random_event((1-ρ)*ζ) #probablity = probability of not going to the ICU * infection probability
        if event == 1
            curr = p[2]
            delete!(E, curr)
            push!(I, curr)
        end
    end
=#
    for p in enumerate(L) #L-->ICU (ICU admission)
        event = generate_random_event(η)
        if event == 1
            curr = p[2]
            delete!(Lcopy, curr)
            push!(ICUcopy, curr)
        end
    end

    for p in enumerate(ICU) #ICU-->R (ICU recovery/death)
        event = generate_random_event(κ)
        if event == 1
            curr = p[2]
            delete!(ICUcopy, curr)
            push!(Rcopy, curr)
        end
    end

    for p in enumerate(I) #I-->R (recovery/death)
        event = generate_random_event(γ)
        if event == 1 #recover/death
            curr = p[2]
            delete!(Icopy, curr)
            push!(Rcopy, curr)
        end
    end

    @assert length(Scopy) + length(Ecopy) + length(Lcopy) + length(ICUcopy) + length(Icopy) + length(Rcopy) == N
    #=
    println("N=", N,", S=", length(S),", E=", length(E), ", L=", length(L),
    ", ICU=", length(ICU), ", I=", length(I),", R=",length(R))
    =#
    return Scopy, Ecopy, Lcopy, ICUcopy, Icopy, Rcopy
end

function move(S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals, S, E, L, ICU, I, R, N, r0, ρ, ζ, η, κ, γ)
    #=
    Function for all sets to be updated and moved for 1 entire cycle in the parameter space
    @params S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals: output vectors to store set sizes
    @params S, E, L, ICU, I, R, N: sets to be moved
    =#
    for t in 1:time_steps #1 period/round of movement e.g. 1 day, 1 month etc.
        S, E, L, ICU, I, R = infect(S, E, L, ICU, I, R, N, r0, ρ, ζ, η, κ, γ) #update sets (e.g. exposure, infection etc.)

        S = move_set(S)
        E = move_set(E)
        L = move_set(L)
        ICU = move_set_ICU(ICU) #ICU movement is different
        I = move_set(I)
        R = move_set(R)

        #update output vectors
        push!(S_vals, length(S)); push!(I_vals, length(I)); push!(R_vals, length(R))
        push!(E_vals, length(E)); push!(L_vals, length(L)); push!(ICU_vals, length(ICU))

    end
    return S_vals, I_vals, R_vals, E_vals, L_vals, ICU_vals, S, E, L, ICU, I, R
end

function move_set(set)
    #=
    Function to move individuals in a set (S, I, R etc.) in the parameter space following the circular paths
    @param set: the set to be moved
    @return: the set with new locations stored
    =#
    set_new = deepcopy(set) #create an independent copy of the set
    for i in enumerate(set) # for every individual
        p = i[2]
        a = p[7]; b = p[8]; r = p[6]; x0 = p[1]; y0 = p[2] #get the parameters

        x_prime = x0 - a; y_prime = y0 - b #treat the center of the circle as the origin for now (i.e. move the origin to the center of the circle, will convert back later)
        xy = SVector(x_prime, y_prime) #create a vector
        rθ = p_from_c(xy) #convert xy to polar
        @assert rθ.r ≈ r #double check conversion
        θ0 = rθ.θ

        rθ_new = Polar(r, θ0 + delta) #update polar angle with the increment value for circluar movement
        xy_new = c_from_p(rθ_new) #convert polar (r,θ0 + delta) to Cartesian (x, y), i.e. the new (x, y) after 1 circular movement step
        x = xy_new[1] + a; y = xy_new[2] + b #back conversion using (a,b), see line 288

        #update x, y in set
        delete!(set_new, p)
        p[1] = x; p[2] = y #set underlying location
        if length(p)>8; pop!(p); pop!(p); end #FIXME: remove previous visualized location
        #FIXME: consider boundary cases and (underlying location vs visualized location)
        #Boundary cases are handled here
        if x < 0.0; x = 0.0
        elseif x > loc_l; x = loc_l; end
        if y < 0; y = 0.0
        elseif y > loc_l; y = loc_l; end

        push!(p, x); push!(p, y) #add visualized location (for boundary cases)
        push!(set_new, p) # add the updated person array back to the set
    end
    return set_new
end

function move_set_ICU(ICU)
    #=
    Function to have individuals in the ICU group stay at the hospital in the parameter space
    Note that the underlying locations of the ICU patients are not changed, which means they would return to the admitted location after recovery
    @param set: the set to be moved
    @return: the set with new locations stored
    =#
    set_new = deepcopy(ICU)
    for i in enumerate(ICU) #for every person in the ICU
        p = i[2]
        a = p[7]; b = p[8]; r = p[6]; x0 = p[1]; y0 = p[2] #get parameters
        delete!(set_new, p)
        if length(p)>8
            p[9]= hospital_x; p[10]=hospital_y #update visualized location
        else
            push!(p, hospital_x); push!(p, hospital_y) #update visualized location
        end
        push!(set_new, p) # add the updated person array back to the set
    end
    return set_new
end
# return the standard deviation
function stdDev(p, n)
    return (p * (1-p) / n)^.5
end
