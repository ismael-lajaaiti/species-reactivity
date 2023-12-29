module SpeciesReactivity

using DifferentialEquations
using Distributions
using LinearAlgebra
using QuadGK

include("community.jl")
include("perturbation.jl")
include("reactivity.jl")

export Community
export amplitude
export assemble!
export complexity
export create_communities
export departure_from_normality
export equilibrium_abundance
export expected_proportional_intensity
export expected_reactivity_squared
export get_reactivity
export is_stable
export isotrope_perturbation
export jacobian
export jacobian_dynamics
export keep_surviving!
export lotka_volterra
export nonlinear_excursion_data
export nonlinear_reactivity
export proportional_perturbation
export random_competition_matrix
export remove_diagonal
export response
export richness
export shannon_diversity
export trajectory_error
export transitional_regime_duration

end
