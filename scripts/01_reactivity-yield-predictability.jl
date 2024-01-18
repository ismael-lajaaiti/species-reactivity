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

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n_communities = 1 # Number of communities.
mu = 0.0
sigma = 0.1
t_max = 1_000 # Duration of observation.
create_interaction_matrix(S) = random_competition_matrix(S, mu, sigma)
communities = create_communities(S, n_communities; create_interaction_matrix)[S]
create_df() = DataFrame(;
    community_id = Int64[],
    yield = Float64[],
    reactivity = Float64[],
    overall_response = Float64[],
    predictability = Float64[],
)
n_perturbations = 1_000
df_to_merge = Any[nothing for _ in 1:n_communities]
Threads.@threads for community_idx in 1:n_communities
    com = communities[community_idx]
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
            push!(
                df,
                [community_idx, yield[i], reactivity[i], overall_response, predictability],
            )
        end
    end
    df_to_merge[community_idx] = df
end
df = reduce(vcat, df_to_merge)

processed_df = combine(
    groupby(df, [:reactivity, :yield, :community_id]),
    :predictability => mean => :predictability,
    :predictability => (x -> std(x) / sqrt(n_perturbations)) => :predictability_error,
    :overall_response => mean => :overall_response,
    :overall_response =>
        (x -> std(x) / sqrt(n_perturbations)) => :overall_response_error,
)

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
    for (i, sdf) in enumerate(groupby(processed_df, :community_id))
        com = communities[i]
        eta_eq = equilibrium_abundance(com)
        reactivity = [get_reactivity(com.A, eta_eq, i) for i in 1:S]
        exp_reactivity_full = sqrt.(expected_reactivity_squared.(eta, Ref(com)))
        exp_reactivity_naive = sqrt.(expected_reactivity_squared_naive.(eta, Ref(com)))
        scatter!(eta_eq, reactivity; color, markersize, strokewidth)
        linewidth = 2
        alpha = 0.7
        naive_exp = lines!(eta, exp_reactivity_naive; label = "Naive", linewidth, alpha)
        full_exp = lines!(eta, exp_reactivity_full; label = "Full", linewidth, alpha)
    end
    axislegend()
    ax2 = Axis(
        fig[1, 2];
        xlabel = "Species reactivity",
        ylabel = "Species overall \n response intensity",
        yscale = log10,
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
    save(
        # "figures/01_reactivity-yield-overall-response.pdf",
        "/tmp/plot.png",
        fig;
        size = (width, height),
        pt_per_unit = 1,
    )
end

# Supplementary figure.
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
    save(
        "figures/SI_01_reactivity-yield-predictability.pdf",
        # "/tmp/plot.png",
        fig;
        size = (width, height),
        pt_per_unit = 1,
    )
end
