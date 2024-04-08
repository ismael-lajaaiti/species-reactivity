using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using LinearAlgebra
using Graphs
using QuadGK
using Random
using SpeciesReactivity
using Statistics

Random.seed!(111) # For reproduciblity.

include("makie-theme.jl")

function random_matrix(S, mu, sigma)
    A = rand(Normal(mu, sigma), S, S)
    A - Diagonal(A) - I
end

function competition_from_adjacency(A)
    S = size(A, 1)
    matrix = zeros(Float64, S, S)
    for i in 1:S, j in 1:S
        matrix[i, j] = A[i, j] == 1 ? -abs(rand(Normal(0, 0.1))) : 0
    end
    matrix - Diagonal(matrix) - I
end
function competition_er(S)
    g = erdos_renyi(S, 0.6)
    A = adjacency_matrix(g)
    competition_from_adjacency(A)
end
function competition_sbm(S)
    g = stochastic_block_model([24 8; 8 24], [ceil(Int, S / 2), floor(Int, S / 2)])
    A = adjacency_matrix(g)
    competition_from_adjacency(A)
end

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
sigma = 0.1
t_max = 1_000 # Duration of observation.
n_perturbations = 100
create_df() = DataFrame(;
    yield = Float64[],
    reactivity = Float64[],
    overall_response = Float64[],
    g_type = String[],
)
com_dict = Dict()
df = create_df()
network_dict = Dict("Erdős–Rényi" => competition_er, "SBM" => competition_sbm)
for (network_type, create_interaction_matrix) in network_dict
    @info "Graph: $network_type"
    com = create_communities(S, 1; create_interaction_matrix)[S][1]
    intensity = 0.6
    A = com.A
    com_dict[network_type] = com
    yield = equilibrium_abundance(com)
    reactivity = [get_reactivity(com.A, yield, i) for i in 1:S]
    for k in 1:n_perturbations
        @info "Perturbation $k"
        x0 = proportional_perturbation(yield, intensity * sqrt(S), 1, true)
        r = SpeciesReactivity.response(com, x0)
        for i in 1:S
            overall_response =
                quadgk(t -> abs(r.nonlinear(t; idxs = i)) / yield[i], 0, t_max)[1]
            push!(df, [yield[i], reactivity[i], overall_response, network_type])
        end
    end
end

processed_df = combine(
    groupby(df, [:reactivity, :yield, :g_type]),
    :overall_response => mean => :overall_response,
    # :overall_response => maximum => :overall_response_worst,
    :overall_response =>
        (x -> std(x) / sqrt(n_perturbations)) => :overall_response_error,
)

with_theme(publication_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    c = fig[2, 1] = GridLayout()
    d = fig[2, 2] = GridLayout()
    e = fig[3, 1] = GridLayout()
    f = fig[3, 2] = GridLayout()
    markersize = 8
    strokewidth = 0.5
    color = :grey
    for (i, g_type) in enumerate(keys(network_dict))
        sub_df = filter(row -> row.g_type == g_type, processed_df)
        com = com_dict[g_type]
        A = com.A
        ax = Axis(
            fig[3, i];
            xlabel = "Species reactivity",
            ylabel = i == 1 ? "Species response \n intensity" : "",
            yscale = log10,
        )
        errorbars!(
            sub_df.reactivity,
            sub_df.overall_response,
            sub_df.overall_response_error;
            linewidth = 1,
            whiskerwidth = 3,
            color = :black,
        )
        scatter!(sub_df.reactivity, sub_df.overall_response; color, markersize, strokewidth)
        ax = Axis(
            fig[2, i];
            xlabel = "Species relative yield",
            ylabel = i == 1 ? "Species reactivity" : "",
        )
        eta_min = minimum(sub_df.yield)
        eta_max = maximum(sub_df.yield)
        eta = range(eta_min, eta_max, 100)
        A_nodiag = A - Diagonal(A)
        C = round(count(!=(0), A_nodiag) / (S * (S - 1)); digits = 2)
        exp_reactivity_full = sqrt.(expected_reactivity_squared.(eta, Ref(com)))
        scatter!(sub_df.yield, sub_df.reactivity; color, markersize, strokewidth)
        lines!(
            eta,
            exp_reactivity_full;
            label = "Expectation",
            linewidth = 2,
            alpha = 0.7,
            color = palette[5],
        )
        ax = Axis(fig[1, i]; title = g_type * " (C = $C)")
        heatmap!(A_nodiag; colormap = :greys)
    end
    for (layout, label) in zip([a, b, c, d, e, f], ["A", "B", "C", "D", "E", "F"])
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
    height = 3 * width * 0.7 / width_height_ratio
    save_figure(
        "figures/SI_network-structure",
        # "/tmp/plot",
        fig,
        (width, height),
    )
end

