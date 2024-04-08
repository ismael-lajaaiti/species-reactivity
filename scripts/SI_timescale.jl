using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using LinearAlgebra
using QuadGK
using Random
using SpeciesReactivity
using Statistics

Random.seed!(111) # For reproduciblity.

include("makie-theme.jl")

function overall_response(solution, i, yield_i)
    quadgk(t -> abs(solution.nonlinear(t; idxs = i)) / yield_i, 0, 1_000)[1] / 1_000
end
function overall_response(solution, yield)
    [overall_response(solution, i, yield_i) for (i, yield_i) in enumerate(yield)]
end

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
n_rep = 1_000
S = 50 # Community richness.
n_communities = 1 # Number of communities.
mu = 0.0
sigma = 0.1
create_interaction_matrix(S) = random_competition_matrix(S, mu, sigma)
create_growth_rates(S) = rand(Uniform(0.1, 1), S)
c =
    create_communities(S, n_communities; create_interaction_matrix, create_growth_rates)[S][1]
# c.r[1] = 10
yield = equilibrium_abundance(c)
reactivity = [get_reactivity(c.A, yield, i) for i in 1:S]
df = DataFrame(; reactivity = Float64[], overall_response = Float64[])
for k in 1:n_rep
    @info "Repetition $k"
    x0 = proportional_perturbation(equilibrium_abundance(c), 0.6 * sqrt(S), 1, true)
    solution = response(c, x0)
    if all(abs.(solution.nonlinear(10_000)) .< 1e-6)
        for (i, yield_i) in enumerate(yield)
            R_sp = overall_response(solution, i, yield_i) # / abs(x0[i])
            push!(df, (reactivity[i], R_sp))
        end
    end
end
error(x) = std(x) / sqrt(length(x))
df_mean = combine(groupby(df, :reactivity), :overall_response => mean, :overall_response => error)

with_theme(publication_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    ax1 = Axis(
        fig[1, 1];
        xlabel = "Species reactivity",
        ylabel = "Species response intensity",
        yscale = log10,
    )
    errorbars!(
        c.r .* df_mean.reactivity,
        df_mean.overall_response_mean,
        df_mean.overall_response_error;
        linewidth = 1,
        whiskerwidth = 3,
        color = :black,
    )
    scatter!(c.r .* reactivity, df_mean.overall_response_mean, color = :grey, strokewidth = 1)
    ax2 = Axis(
        fig[1, 2];
        xlabel = "Species reactivity",
        ylabel = "",
        yscale = log10,
    )
    errorbars!(
        df_mean.reactivity,
        c.r .* df_mean.overall_response_mean,
        df_mean.overall_response_error;
        linewidth = 1,
        whiskerwidth = 3,
        color = :black,
    )
    scatter!(reactivity, c.r .* df_mean.overall_response_mean, color = :grey, strokewidth = 1)
    isdir("figures") || mkdir("figures")
    width = full_page_width * cm_to_pt
    height = width * 0.7 / width_height_ratio
    save(
        "figures/SI_01_timescale.pdf",
        # "/tmp/plot.png",
        fig;
        size = (width, height),
        pt_per_unit = 1,
    )
end
