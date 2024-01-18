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

include("makie-theme.jl")

Random.seed!(113) # For reproducibility.

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
        if richness(c) >= 6
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
                cond_pos = (se == :positive && p < 0.05 && cor(K, eta) > 0)
                cond_neg = (se == :negative && p < 0.05 && cor(K, eta) < 0)
                if cond_null || cond_pos || cond_neg
                    push!(communities, c)
                    push!(K_communities, K)
                    push!(eta_communities, eta)
                    iter += 1
                end
            end
        end
    end
    communities, K_communities, eta_communities
end

function get_K(c, surviving_species, df_K)
    df = df_K[in.(df_K[:, "Species"], Ref(surviving_species)), :]
    [rand(Normal(row["50"], row["sigma"])) for row in eachrow(df)]
end

communities, K_communities, eta_communities =
    generate_communities_from_data(10, df_A, df_K; se = :any)
eta_flat = reduce(vcat, eta_communities)
eta_linrange = LinRange(0, 0.7, 100)
reactivity_communities = get_reactivity.(communities)
abundance_communities = [eta .* K for (eta, K) in zip(eta_communities, K_communities)]
abundance_communities = [abundance ./ sum(abundance) for abundance in abundance_communities]
se_communities = [cov(K, eta) for (K, eta) in zip(K_communities, eta_communities)]

n_perturbations = 100
n_communities = 100
intensity = 0.6
t_max = 1_000
se_dict_value = Dict(:negative => -1, :null => 0, :positive => 1)
se_symbols = [:negative, :null, :positive]
create_df() = DataFrame(; selection_effect = Int64[], community_response = Float64[])
df_vec = Any[nothing, nothing, nothing]
Threads.@threads for se_idx in 1:3
    se = se_symbols[se_idx]
    df = create_df()
    coms, K_coms, eta_coms = generate_communities_from_data(n_communities, df_A, df_K; se)
    for (com, K, eta) in zip(coms, K_coms, eta_coms)
        @info "New community."
        S = richness(com)
        A = com.A
        for k in 1:n_perturbations
            x0 = proportional_perturbation(eta, intensity * sqrt(S), 1, true)
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
with_theme(p_theme) do
    fig = Figure()
    panel_a = fig[1, 1] = GridLayout()
    panel_b = fig[2, 1] = GridLayout()
    panel_c = fig[3, 1] = GridLayout()
    color = :grey
    markersize = 9
    linewidth = 2
    strokewidth = 0.5
    # Panel A: reactivity vs. relative yield.
    ax1 = Axis(fig[1, 1]; xlabel = "Species relative yield", ylabel = "Species reactivity")
    for (eta, reactivity) in zip(eta_communities, reactivity_communities)
        scatter!(eta, reactivity; markersize, strokewidth)
    end
    for (eta, com) in zip(eta_communities, communities)
        eta_linrange = LinRange(minimum(eta), maximum(eta), 100)
        exp_reactivity = sqrt.(expected_reactivity_squared.(eta_linrange, Ref(com)))
        lines!(eta_linrange, exp_reactivity; linewidth, alpha = 0.6)
    end
    # Panel B: reactivity vs. abundance.
    ax2 = Axis(fig[2, 1]; xlabel = "Species abundance", ylabel = "Species reactivity")
    for (N, reactivity) in zip(abundance_communities, reactivity_communities)
        scatter!(N, reactivity; markersize, strokewidth)
    end
    reactivity_flat = reduce(vcat, reactivity_communities)
    abundance_flat = reduce(vcat, abundance_communities)
    df_tmp = DataFrame(; reactivity_flat, abundance_flat)
    ols = lm(@formula(reactivity_flat ~ abundance_flat), df_tmp)
    p_value = ftest(ols.model).pval
    a, b = coef(ols)
    xlims = [minimum(abundance_flat), maximum(abundance_flat)]
    ylims = a .+ b .* xlims
    linestyle = p_value > 0.05 ? :dot : :solid
    lines!(xlims, ylims; color = :grey, linestyle, linewidth = 2)
    # Panel C: community response vs. selection effect.
    ax3 = Axis(
        fig[3, 1];
        xlabel = "Selection effect",
        ylabel = "Log community overall \n response intensity",
        xticks = (-1:1, ["Negative", "Null", "Positive"]),
    )
    for (se, color) in zip([-1, 0, 1], palette[[1, 3, 2]])
        sdf = df[df.selection_effect.==se, :]
        boxplot!(
            sdf.selection_effect,
            log10.(sdf.community_response);
            show_outliers = false,
            whiskerwidth = 0.5,
            color,
        )
    end
    # Label panels.
    for (layout, label) in zip([panel_a, panel_b, panel_c], ["A", "B", "C"])
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
    height = width * 3 / width_height_ratio
    save(
        # "/tmp/plot.png",
        "figures/03_data.pdf",
        fig;
        size = (width, height),
        pt_per_unit = 1,
    )
end

# Figure for supplementary material.
with_theme(p_theme) do
    fig = Figure()
    markersize = 9
    linewidth = 2
    strokewidth = 0.5
    panel_a = fig[1, 1] = GridLayout()
    panel_b = fig[1, 2] = GridLayout()
    # Panel A: Abundance vs. relative yield.
    ax = Axis(
        fig[1, 1];
        xlabel = "Relative yield",
        ylabel = "Abundance",
        title = L"\mathrm{Cov}(N, \eta) > 0",
    )
    for (eta, N) in zip(eta_communities, abundance_communities)
        scatter!(eta, N; markersize, strokewidth)
    end
    # Panel B: carrying capacity vs. relative yield (selection effect).
    ax1 = Axis(
        fig[1, 2];
        xlabel = "Relative yield",
        ylabel = "Carrying capacity",
        title = "Overall \n selection effet ≃ 0",
    )
    for (eta, K) in zip(eta_communities, K_communities)
        scatter!(eta, K; markersize, strokewidth)
    end
    K_flat = reduce(vcat, K_communities)
    df_tmp = DataFrame(; eta_flat, K_flat)
    ols = lm(@formula(K_flat ~ eta_flat), df_tmp)
    p_value = ftest(ols.model).pval
    a, b = coef(ols)
    xlims = [minimum(eta_flat), maximum(eta_flat)]
    ylims = a .+ b .* xlims
    linestyle = p_value > 0.05 ? :dot : :solid
    lines!(xlims, ylims; color = :grey, linestyle, linewidth = 3)
    # Label panels.
    for (layout, label) in zip([panel_a, panel_b], ["A", "B"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    isdir("figures") || mkdir("figures")
    save(
        "/tmp/plot.png",
        # "figures/SI_03_data.png",
        fig;
        size = (200 * 1.3 * 2, 200),
        px_per_unit = 3,
    )
end
