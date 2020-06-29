using Distributions
function sir_step(S,I,R, β, dt, )
end

function sir_init(N, eta)
    S = round(N*eta)
    I = 1
    R = round(N*(1-eta))
end
    
