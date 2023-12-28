module RareLotkaVolterra

using CairoMakie
using DifferentialEquations
using Distributions
using LinearAlgebra
import Makie: texfont
using QuadGK

include("community.jl")
include("perturbation.jl")
include("makie-theme.jl")

export Community
export amplitude
export assemble!
export complexity
export create_communities
export departure_from_normality
export equilibrium_abundance
export expected_proportional_intensity
export is_stable
export isotrope_perturbation
export jacobian
export jacobian_dynamics
export keep_surviving!
export lotka_volterra
export nonlinear_excursion_data
export nonlinear_reactivity
export proportional_perturbation
export publication_theme
export random_competition_matrix
export response
export richness
export shannon_diversity
export trajectory_error
export transitional_regime_duration

end
