using CairoMakie
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics

# Create communities.
S = 20 # Species richness.
sigma = 0.1 # Interaction strength.
n = 10 # Number of communities.
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
communities = create_communities(S, n; create_interaction_matrix, verbose)[S]

# Create perturbations.
n_perturbation = 100
function isotrope_perturbation()
    x = rand(Normal(), S)
    x / norm(x)
end
perturbation_vec = [isotrope_perturbation() for _ in 1:n_perturbation]
intensity_values = [0.1, 0.2, 0.5]

for com in communities, x0 in perturbation_vec
    Neq = equilibrium_abundance(com)
    Ntot = sum(Neq)
    for i0 in intensity_values
        no_extinction = all(Neq + i0 * x0 .> 0) # True if there is no extinction.
    end
end
