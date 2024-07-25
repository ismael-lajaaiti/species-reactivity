import SpecialFunctions: erfinv
using CSV
using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using GLM
using HypothesisTests
using LinearAlgebra
using QuadGK
using Random
using SpeciesReactivity
using Statistics
using StatsBase
include("makie-theme.jl")

Random.seed!(112) # For reproducibility.

# Read data.
df_K = CSV.read("../data/carrying_capacity.csv", DataFrame)
df_A = CSV.read("../data/interactions.csv", DataFrame)

infer_std(p, q, mu) = (q - mu) / (sqrt(2) * erfinv(2p - 1))

# Process interaction dataframe.
transform!(
    df_A,
    "Species i" => ByRow(x -> x + 1),
    "Species j" => ByRow(x -> x + 1);
    renamecols = false,
)
transform!(df_A, ["10", "50"] => ByRow((q, mu) -> infer_std(0.1, q, mu)) => "sigma10")
transform!(df_A, ["90", "50"] => ByRow((q, mu) -> infer_std(0.9, q, mu)) => "sigma90")
transform!(df_A, ["sigma10", "sigma90"] => ByRow((x, y) -> mean([x, y])) => "sigma")

# Process carrying capacity dataframe.
transform!(df_K, "Species" => ByRow(x -> x += 1); renamecols = false)
transform!(df_K, ["10", "50"] => ByRow((q, mu) -> infer_std(0.1, q, mu)) => "sigma10")
transform!(df_K, ["90", "50"] => ByRow((q, mu) -> infer_std(0.9, q, mu)) => "sigma90")
transform!(df_K, ["sigma10", "sigma90"] => ByRow((x, y) -> mean([x, y])) => "sigma")

function create_interaction_matrix(df_A)
    @assert maximum(df_A[:, "Species i"]) == maximum(df_A[:, "Species j"])
    S = maximum(df_A[:, "Species i"])
    A = zeros(S, S) - I
    for row in eachrow(df_A)
        i = row["Species i"]
        j = row["Species j"]
        mu = row["50"]
        sigma = row["sigma"]
        i != j && (A[i, j] = -rand(Normal(mu, sigma)))
    end
    A
end

function generate_communities_from_data(n_communities, df_A, df_K; se = :any)
    r = fill(1, 8)
    iter = 0
    communities = []
    eta_communities = []
    K_communities = []
    while length(communities) < n_communities || iter > 10_000
        A = create_interaction_matrix(df_A)
        c = Community(A, r)
        idx = assemble!(c; return_surviving = true)
        if richness(c) >= 6 && all(equilibrium_abundance(c) .> 0)
            eta = equilibrium_abundance(c)
            K = get_K(c, idx, df_K)
            if se == :any
                push!(communities, c)
                push!(K_communities, K)
                push!(eta_communities, eta)
                iter += 1
            else
                p = pvalue(CorrelationTest(K, eta))
                cond_null = (se == :null && p > 0.05)
                cond_pos = (se == :positive && p < 0.01 && cor(K, eta) > 0)
                cond_neg = (se == :negative && p < 0.01 && cor(K, eta) < 0)
                if cond_null || cond_pos || cond_neg
                    push!(communities, c)
                    push!(K_communities, K)
                    push!(eta_communities, eta)
                    iter += 1
                    @info length(communities)
                end
            end
        end
    end
    communities, K_communities, eta_communities
end

function get_K(c, surviving_species, df_K)
    df = df_K[in.(df_K[:, "Species"], Ref(surviving_species)), :]
    K = fill(0, richness(c))
    while any(K .<= 0)
        K = [rand(Normal(row["50"], row["sigma"])) for row in eachrow(df)]
    end
    K
end

Random.seed!(123) # For reproducibility.
communities, K_communities, eta_communities =
    generate_communities_from_data(10, df_A, df_K; se = :any)
eta_flat = reduce(vcat, eta_communities)
eta_linrange = LinRange(0, 0.7, 100)
reactivity_communities = get_reactivity.(communities)
abundance_communities = [eta .* K for (eta, K) in zip(eta_communities, K_communities)]
abundance_communities = [abundance ./ sum(abundance) for abundance in abundance_communities]
se_communities = [cov(K, eta) for (K, eta) in zip(K_communities, eta_communities)]

# Generate communities in parallel for time efficiency.
n_threads = 5
n_com_per_thread = 10
community_set2 = Dict()
for se in se_symbols
    @info se
    c_list = Any[undef for _ in 1:n_threads]
    K_list = Any[undef for _ in 1:n_threads]
    eta_list = Any[undef for _ in 1:n_threads]
    Threads.@threads for thr in 1:n_threads
        c, K, eta = generate_communities_from_data(n_com_per_thread, df_A, df_K; se)
        c_list[thr] = c
        K_list[thr] = K
        eta_list[thr] = eta
    end
    community_set2[se] = (vcat(c_list...), vcat(K_list...), vcat(eta_list...))
end

n_perturbations = 100
n_communities = 10
intensity = 0.9
t_max = 1_000
se_dict_value = Dict(:negative => -1, :null => 0, :positive => 1)
se_symbols = [:positive, :negative]
community_set = Dict(
    se => generate_communities_from_data(n_communities, df_A, df_K; se) for
    se in se_symbols
)

create_df() = DataFrame(; selection_effect = Int64[], community_response = Float64[])
df_vec = Any[nothing for _ in 1:length(se_symbols)]
Threads.@threads for se_idx in 1:2
    se = se_symbols[se_idx]
    df = create_df()
    coms, K_coms, eta_coms = community_set2[se]
    for (com, K, eta) in zip(coms, K_coms, eta_coms)
        @info "New community."
        S = richness(com)
        A = com.A
        for k in 1:n_perturbations
            x0 = prop_perturbation(eta, intensity; no_extinction = true)
            # x0 = proportional_perturbation(eta, intensity * sqrt(S), 1, true)
            r = SpeciesReactivity.response(com, x0)
            delta_total_abundance(delta_eta) = abs(sum(K .* delta_eta))
            a = delta_total_abundance(x0) * t_max
            com_response =
                quadgk(t -> (delta_total_abundance(r.nonlinear(t)) / a), 0, t_max)[1]
            push!(df, [se_dict_value[se], com_response])
        end
    end
    df_vec[se_idx] = df
end
df = reduce(vcat, df_vec)
processed_df = combine(
    groupby(df, :selection_effect),
    :community_response => median => :community_response,
    :community_response =>
        (x -> std(x) / sqrt(n_perturbations)) => :community_response_err,
)

# Figure for main manuscript.
with_theme(publication_theme) do
    fig = Figure()
    panel_a = fig[1, 1:2] = GridLayout()
    panel_b = fig[2, 1] = GridLayout()
    panel_c = fig[2, 2] = GridLayout()
    panel_d = fig[3, 1:2] = GridLayout()
    color = :grey
    markersize = 9
    linewidth = 2
    strokewidth = 0.5
    # Panel A: reactivity vs. relative yield.
    ax1 =
        Axis(fig[1, 1:2]; xlabel = "Species relative yield", ylabel = "Species reactivity")
    for (eta, reactivity) in zip(eta_communities, reactivity_communities)
        scatter!(eta, reactivity; markersize, strokewidth)
    end
    for (eta, com) in zip(eta_communities, communities)
        eta_linrange = LinRange(minimum(eta), maximum(eta), 100)
        exp_reactivity = sqrt.(expected_reactivity_squared.(eta_linrange, Ref(com)))
        lines!(eta_linrange, exp_reactivity; linewidth, alpha = 0.6)
    end
    # Panel B: reactivity vs. abundance.
    title_se = ["Positive \n selection effect", "Negative \n selection effect"]
    for (i, se) in enumerate(se_symbols)
        ylabel = i == 1 ? "Species reactivity" : ""
        ax = Axis(
            fig[2, i];
            xlabel = "Species abundance",
            ylabel,
            # yscale = log10,
            xscale = log10,
            title = title_se[i],
            titlesize = 12,
        )
        R0_N_slopes = Float64[]
        com_set, K_set, eta_set = community_set2[se]
        reactivity_set = get_reactivity.(com_set)
        abundance_set = [eta .* K for (eta, K) in zip(eta_set, K_set)]
        abundance_set = [abundance ./ sum(abundance) for abundance in abundance_set]
        rho = round(cor(vcat(reactivity_set...), vcat(abundance_set...)); digits = 2)
        text!(
            ax,
            0.05,
            0.85;
            text = "ρ = $rho",
            offset = (3, -2),
            align = (:left, :top),
            fontsize = 10,
        )
        for (abundance, reactivity, K) in zip(abundance_set, reactivity_set, K_set)
            scatter!(abundance, reactivity; markersize = 7, strokewidth, alpha = 1)
        end
        i == 1 && (ax.xticks = [1e-2, 1e-1, 1])
        i == 2 && (ax.xticks = [1e-2, 0.05, 0.3])
    end
    # Panel C: community response vs. selection effect.
    ax3 = Axis(
        fig[3, 1:2];
        xlabel = "Log community overall \n response intensity",
        ylabel = "Density",
    )
    for (se, color) in zip([1, -1], palette[[2, 1]])
        sdf = df[df.selection_effect.==se, :]
        tmp_hist = density!(
            ax3,
            log10.(sdf.community_response);
            color = (color, 0.3),
            strokecolor = color,
            strokewidth = 1.4,
            label = se == 1 ? "Positive" : "Negative",
        )
    end
    axislegend("Selection effect"; position = :rt, rowgap = 0)
    # Label panels.
    for (layout, label) in zip([panel_a, panel_b, panel_c, panel_d], letters)
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    isdir("figures") || mkdir("figures")
    width = two_third_page_width * cm_to_pt
    height = width * 2.6 / width_height_ratio
    save_figure(
        # "/tmp/plot",
        "figures/03_data",
        fig,
        (width, height),
    )
end
