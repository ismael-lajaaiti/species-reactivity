using RareLotkaVolterra

S = 10
community = create_communities(S, 1)[S][1]
Neq = equilibrium_abundance(community)
rarest = argmin(Neq)
richest = argmax(Neq)
i0 = 0.1
x0 = zeros(Float64, S)
x0[rarest] = i0
r = response(community, x0)
