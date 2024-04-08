using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using LinearAlgebra
using QuadGK
using Random
using SpeciesReactivity
using Statistics

Random.seed!(123) # For reproduciblity.

include("makie-theme.jl")

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n_communities = 1 # Number of communities.
mu = 0.0
sigma = 0.1
t_max = 1_000 # Duration of observation.
create_interaction_matrix(S) = random_competition_matrix(S, mu, sigma)
com = create_communities(S, n_communities; create_interaction_matrix)[S][1]
create_df() = DataFrame(; community_id = Int64[], yield = Float64[], reactivity = Float64[])
n_perturbations = 100 # For SI Fig. set to 1_000.
intensity = 0.6
df_to_merge = Any[nothing for _ in 1:n_communities]
df = create_df()
A = com.A
A_nodiag = A - Diagonal(A)
yield = equilibrium_abundance(com)
reactivity = [get_reactivity(A, yield, i) for i in 1:S]
mean_a_ij_square = vec(mean(A .^ 2; dims = 2))


# Main figure.
eta = LinRange(0.0, 0.5, 100)
with_theme(publication_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    markersize = 8
    strokewidth = 0.5
    color = :grey
    ax1 = Axis(fig[1, 1]; xlabel = "Species relative yield", ylabel = "Species reactivity")
    eta_eq = equilibrium_abundance(com)
    reactivity = [get_reactivity(com.A, eta_eq, i) for i in 1:S]
    exp_reactivity_full = sqrt.(expected_reactivity_squared.(eta, Ref(com)))
    exp_reactivity_naive = sqrt.(expected_reactivity_squared_naive.(eta, Ref(com)))
    scatter!(eta_eq, reactivity; color, markersize, strokewidth)
    linewidth = 2
    alpha = 0.7
    naive_exp = lines!(
        eta,
        exp_reactivity_naive;
        label = "Naive",
        linewidth,
        alpha,
        color = palette[6],
    )
    full_exp = lines!(
        eta,
        exp_reactivity_full;
        label = "Full",
        linewidth,
        alpha,
        color = palette[5],
    )
    axislegend()
    ax2 = Axis(
        fig[1, 2];
        xlabel = "Species relative yield",
        ylabel = "Mean received interaction",
    )
    scatter!(yield, mean_a_ij_square; color, markersize, strokewidth)
    for (layout, label) in zip([a, b], ["A", "B"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    isdir("figures") || mkdir("figures")
    width = full_page_width * cm_to_pt
    height = width * 0.7 / width_height_ratio
    save_figure(
        "figures/SI_expectation",
        # "/tmp/plot.png",
        fig,
        (width, height),
    )
end
