function updateSIR(population)
    S = population[1]
    I = population[2]
    R = population[3]
    V = population[4]

    S_new = S - λ * S * I * dt - p * S *dt
    I_new = I + λ * S * I * dt - γ * I * dt
    R_new = R + γ * I * dt
    V_new = V + p * S * dt
    #p = fraction of Ppl getting successfully vaccinated during time unit

    return [S_new I_new R_new V_new]
end
