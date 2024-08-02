using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using LinearAlgebra
using QuadGK
using Random
using SpeciesReactivity
using Statistics
using StatsBase

Random.seed!(111) # For reproduciblity.

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
create_df() = DataFrame(;
    yield = Float64[],
    reactivity = Float64[],
    overall_response = Float64[],
    predictability = Float64[],
)
n_perturbations = 100
intensity = 0.6
df = create_df()
A = com.A
yield = equilibrium_abundance(com)
reactivity = [get_reactivity(com.A, yield, i) for i in 1:S]

for k in 1:n_perturbations
    @info "Perturbation $k"
    x0 = proportional_perturbation(yield, intensity * sqrt(S), 1, true)
    r = SpeciesReactivity.response(com, x0)
    for i in 1:S
        deviation = trajectory_error(
            t -> r.nonlinear(t; idxs = i),
            t -> r.linear(t; idxs = i);
            tspan = (0, t_max),
        )
        predictability = exp(-deviation)
        overall_response =
            quadgk(t -> abs(r.nonlinear(t; idxs = i)) / yield[i], 0, t_max)[1]
        push!(df, [yield[i], reactivity[i], overall_response, predictability])
    end
end

processed_df = combine(
    groupby(df, [:reactivity, :yield]),
    :predictability => mean => :predictability,
    :predictability => (x -> std(x) / sqrt(n_perturbations)) => :predictability_error,
    :overall_response => mean => :overall_response,
    :overall_response => maximum => :overall_response_worst,
    :overall_response =>
        (x -> std(x) / sqrt(n_perturbations)) => :overall_response_error,
)

# First figure.
with_theme(publication_theme) do
    fig = Figure()
    markersize = 8
    strokewidth = 0.5
    color = :grey
    ax = Axis(fig[1, 1]; xlabel = "Species reactivity", ylabel = "Species predictability")
    errorbars!(
        processed_df.reactivity,
        processed_df.predictability,
        processed_df.predictability_error;
        linewidth = 1,
        whiskerwidth = 3,
        color = :black,
    )
    scatter!(
        processed_df.reactivity,
        processed_df.predictability;
        color,
        markersize,
        strokewidth,
    )
    isdir("figures") || mkdir("figures")
    width = two_third_page_width * cm_to_pt
    height = width / width_height_ratio
    save_figure(
        # "figures/SI_predictability",
        "/tmp/plot",
        fig,
        (width, height),
    )
end

# Second supplementary figure.
eta = LinRange(0.0, 0.5, 100)
with_theme(publication_theme) do
    fig = Figure()
    panel_a = fig[1, 1] = GridLayout()
    panel_b = fig[1, 2] = GridLayout()
    panel_c = fig[1, 3] = GridLayout()
    markersize = 8
    strokewidth = 0.5
    color = :grey
    ax2 = Axis(
        fig[1, 2];
        xlabel = "Species average \n initial response",
        # ylabel = "Species overall response intensity",
        yscale = log10,
        # xscale = log10,
    )
    errorbars!(
        -processed_df.yield,
        processed_df.overall_response,
        processed_df.overall_response_error;
        linewidth = 1,
        whiskerwidth = 3,
        color = :black,
    )
    scatter!(
        -processed_df.yield,
        processed_df.overall_response;
        color,
        markersize,
        strokewidth,
    )
    ax1 = Axis(
        fig[1, 1];
        xlabel = "Species reactivity",
        ylabel = "Species overall \n response intensity",
        yscale = log10,
        # xscale = log10,
    )
    errorbars!(
        processed_df.reactivity,
        processed_df.overall_response,
        processed_df.overall_response_error;
        linewidth = 1,
        whiskerwidth = 3,
        color = :black,
    )
    scatter!(
        processed_df.reactivity,
        processed_df.overall_response;
        color,
        markersize,
        strokewidth,
    )
    ax3 = Axis(
        fig[1, 3];
        xlabel = "Species maximum \n displacement",
        # ylabel = "Species response intensity",
        yscale = log10,
        # xscale = log10,
    )
    ratio = processed_df.reactivity ./ processed_df.yield
    f(x) = x / sqrt(1 + x^2)
    errorbars!(
        f.(ratio),
        processed_df.overall_response,
        processed_df.overall_response_error;
        linewidth = 1,
        whiskerwidth = 3,
        color = :black,
    )
    scatter!(f.(ratio), processed_df.overall_response; color, markersize, strokewidth)
    hideydecorations!(ax2)
    hideydecorations!(ax3)
    for (layout, label) in zip([panel_a, panel_b, panel_c], letters)
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
    height = width * 0.6 / width_height_ratio
    save_figure(
        "figures/SI_reactivity-resilience",
        # "/tmp/plot",
        fig,
        1.3 .* (width, height),
    )
end
