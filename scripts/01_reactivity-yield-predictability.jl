using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using LinearAlgebra
using QuadGK
using SpeciesReactivity
using Statistics

include("makie-theme.jl")

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n_communities = 1 # Number of communities.
mu = 0.0
sigma = 0.1
create_interaction_matrix(S) = random_competition_matrix(S, mu, sigma)
communities = create_communities(S, n_communities; create_interaction_matrix)[S]
create_df() = DataFrame(;
    community_id = Int64[],
    intensity = Float64[],
    yield = Float64[],
    reactivity = Float64[],
    predictability = Float64[],
)
n_perturbations = 1_000
intensity_values = [0.4, 0.6, 0.9]
color_scheme = ColorSchemes.seaborn_pastel
intensity_colors =
    Dict(intensity => color for (intensity, color) in zip(intensity_values, color_scheme))

function expected_reactivity2(eta_i, com)
    S = richness(com)
    d = S - 1
    A = remove_diagonal(com.A)
    eta = equilibrium_abundance(com)
    mu = sum(A) / (d * S)
    sigma = sum(A .^ 2) / (d * S)
    D_i = get_Di(eta, eta_i)
    a =
        (1 / d) * (1 - eta_i)^2 / (norm(eta)^2 - eta_i^2) + ((d - 1) * sigma) / d -
        ((D_i - 1) * mu^2) / d
    # a = ((d-1) * sigma)/d - ((D_i - 1) * mu^2) / d
    # a = (1 / d) * (1 - eta_i)^2 / (norm(eta)^2 - eta_i^2) - ((D_i - 1) * mu^2) / d
    # a
    a * (norm(eta)^2 - eta_i^2)
end

df_to_merge = Any[nothing for _ in 1:n_communities]
Threads.@threads for community_idx in 1:n_communities
    com = communities[community_idx]
    df = create_df()
    A = com.A
    yield = equilibrium_abundance(com)
    reactivity = [get_reactivity(com.A, yield, i) for i in 1:S]
    # reactivity = [get_reactivity(A, yield, i) for i in 1:S]
    for intensity in intensity_values, k in 1:n_perturbations
        @info "Perturbation $k"
        x0 = proportional_perturbation(yield, intensity * sqrt(S), 1, true)
        r = SpeciesReactivity.response(com, x0)
        for i in 1:S
            deviation = trajectory_error(
                t -> r.nonlinear(t; idxs = i),
                t -> r.linear(t; idxs = i);
                tspan = (0, 10_000),
            )
            predictability = exp(-deviation)
            push!(df, [community_idx, intensity, yield[i], reactivity[i], predictability])
        end
    end
    df_to_merge[community_idx] = df
end
df = reduce(vcat, df_to_merge)
processed_df = combine(
    groupby(df, [:reactivity, :yield, :intensity, :community_id]),
    :predictability => mean => :predictability,
    :predictability => (x -> std(x) / sqrt(n_perturbations)) => :predictability_error,
)

eta = LinRange(0.0, 0.5, 100)
# exp_reactivity2 = expected_reactivity2.(eta, Ref(communities[1]))
with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    markersize = 8
    strokewidth = 0.5
    ax1 = Axis(fig[1, 1]; xlabel = "Relative yield", ylabel = "Reactivity")
    for (i, sdf) in enumerate(groupby(processed_df, :community_id))
        # color = color_scheme[i]
        com = communities[i]
        eta_eq = equilibrium_abundance(com)
        reactivity = [get_reactivity(com.A, eta_eq, i) for i in 1:S]
        exp_reactivity_full = sqrt.(expected_reactivity_squared.(eta, Ref(com)))
        scatter!(eta_eq, reactivity; color = :grey, markersize, strokewidth)
        # lines!(eta, exp_reactivity2; color = :grey, label = "approx")
        lines!(
            eta,
            exp_reactivity_full;
            color = :coral,
            label = "expectation",
            linewidth = 2,
            alpha = 0.7,
        )
    end
    axislegend()
    ax2 = Axis(fig[1, 2]; xlabel = "Reactivity", ylabel = "Predictability")
    for ((marker_idx, intensity), sdf) in
        pairs(groupby(processed_df, [:community_id, :intensity]))
        marker = :circle
        color = intensity_colors[intensity]
        errorbars!(
            sdf.reactivity,
            sdf.predictability,
            sdf.predictability_error;
            linewidth = 1,
            whiskerwidth = 3,
        )
        scatter!(
            sdf.reactivity,
            sdf.predictability;
            label = L"\sqrt{\langle z_0^2 \rangle} = %$intensity",
            color,
            marker,
            markersize,
            strokewidth,
        )
        if marker_idx == 1 && intensity == intensity_values[end]
            axislegend(; position = :lb)
            ylims!(0.5, 0.75)
        end
    end
    for (layout, label) in zip([a, b], ["A", "B"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    save(
        # "figures/reactivity-yield-predictability.eps",
        "/tmp/plot.png",
        fig;
        save = (600, 320),
        px_per_unit = 3,
    )
end
