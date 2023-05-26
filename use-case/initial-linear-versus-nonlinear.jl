using CairoMakie
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics

# Create communities.
S = 20 # Species richness.
sigma = 0.1 # Interaction strength.
n = 10 # Number of communities.
verbose = false
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
communities = create_communities(S, n; create_interaction_matrix, verbose)[S]

# Create perturbations.
n_perturbation = 100
function perturbation()
    x = rand(Normal(), S)
    x / norm(x)
end
perturbation_vec = [perturbation() for _ in 1:n_perturbation]
intensity_values = [0.01, 0.1]

c = communities[1]
Neq = equilibrium_abundance(c)
i_min = argmin(Neq)
i_max = argmax(Neq)
for x0 in perturbation_vec
    intensity = 0.01
    x0 *= mean(Neq) * intensity
    linear_term = norm(Neq[i_min] * sum(c.A[i_min, :] .* x0))^2
    nonlinear_term = norm(x0[i_min] * sum(c.A[i_min, :] .* x0))^2
end
