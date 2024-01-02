import SpecialFunctions: erfinv
using CSV
using CairoMakie
using DataFrames
using Distributions
using GLM
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

r = fill(1, 8)
iter = 0
n_commmunities = 10
communities = []
surviving_species = []
while length(communities) < n_commmunities || iter > 10_000
    A = create_interaction_matrix(df_A)
    c = Community(A, r)
    idx = assemble!(c; return_surviving = true)
    if richness(c) >= 6
        push!(communities, c)
        push!(surviving_species, idx)
    end
    global iter += 1
end

function get_K(c, surviving_species, df_K)
    df = df_K[in.(df_K[:, "Species"], Ref(surviving_species)), :]
    [rand(Normal(row["50"], row["sigma"])) for row in eachrow(df)]
end

eta_communities = equilibrium_abundance.(communities)
eta_flat = reduce(vcat, eta_communities)
eta_linrange = LinRange(0, 0.7, 100)
reactivity_communities = get_reactivity.(communities)
K_communities = get_K.(communities, surviving_species, Ref(df_K))
abundance_communities = [eta .* K for (eta, K) in zip(eta_communities, K_communities)]
abundance_communities = [abundance ./ sum(abundance) for abundance in abundance_communities]

idx_com_simu = 6 # rand(1:n_commmunities)
com_simu = communities[idx_com_simu]
eta_simu = eta_communities[idx_com_simu]
K_simu = K_communities[idx_com_simu]
abundance_simu = eta_simu .* K_simu
idx = argmin(eta_simu) # Most reactive species.
x0 = remove_diagonal(com_simu.A)[idx, :] .* eta_simu
r = SpeciesReactivity.response(com_simu, x0)
time_steps = collect(0:0.1:50)
eta_true = [abs.(r.nonlinear(t)) for t in time_steps]
eta_lin = [abs.(r.linear(t)) for t in time_steps]
abundance_true = [sum(eta_t .* K_simu) for eta_t in eta_true] / sum(abundance_simu)
abundance_lin = [sum(eta_t .* K_simu) for eta_t in eta_lin] / sum(abundance_simu)

with_theme(p_theme) do
    fig = Figure()
    ax1 = Axis(fig[1, 1]; xlabel = "Relative yield", ylabel = "Reactivity")
    color = :grey
    markersize = 9
    linewidth = 2
    strokewidth = 0.5
    for (eta, reactivity) in zip(eta_communities, reactivity_communities)
        scatter!(eta, reactivity; markersize, strokewidth)
    end
    for (eta, com) in zip(eta_communities, communities)
        eta_linrange = LinRange(minimum(eta), maximum(eta), 100)
        exp_reactivity = sqrt.(expected_reactivity_squared.(eta_linrange, Ref(com)))
        lines!(eta_linrange, exp_reactivity; linewidth, alpha = 0.6)
    end
    ax2 = Axis(
        fig[2, 2];
        xlabel = "Relative yield",
        ylabel = "Carrying capacity",
        title = "Selection effet ≃ 0",
    )
    for (eta, K) in zip(eta_communities, K_communities)
        scatter!(eta, K; markersize, strokewidth)
    end
    K_flat = reduce(vcat, K_communities)
    df = DataFrame(; eta_flat, K_flat)
    ols = lm(@formula(K_flat ~ eta_flat), df)
    p_value = ftest(ols.model).pval
    a, b = coef(ols)
    xlims = [minimum(eta_flat), maximum(eta_flat)]
    ylims = a .+ b .* xlims
    linestyle = p_value > 0.05 ? :dot : :solid
    lines!(xlims, ylims; color = :skyblue, linestyle, linewidth = 3)
    ax3 = Axis(
        fig[2, 1];
        xlabel = "Relative yield",
        ylabel = "Abundance",
        title = L"\mathrm{Cov}(N, \eta) > 0",
    )
    for (eta, N) in zip(eta_communities, abundance_communities)
        scatter!(eta, N; markersize, strokewidth)
    end
    ax4 = Axis(fig[3, 1]; xlabel = "Abundance", ylabel = "Reactivity")
    for (N, reactivity) in zip(abundance_communities, reactivity_communities)
        scatter!(N, reactivity; markersize, strokewidth)
    end
    reactivity_flat = reduce(vcat, reactivity_communities)
    abundance_flat = reduce(vcat, abundance_communities)
    df = DataFrame(; reactivity_flat, abundance_flat)
    ols = lm(@formula(reactivity_flat ~ abundance_flat), df)
    p_value = ftest(ols.model).pval
    a, b = coef(ols)
    xlims = [minimum(abundance_flat), maximum(abundance_flat)]
    ylims = a .+ b .* xlims
    linestyle = p_value > 0.05 ? :dot : :solid
    lines!(xlims, ylims; color = :skyblue, linestyle, linewidth = 3)
    ax4 = Axis(fig[4, 1]; xlabel = "Δ tot. abundance", ylabel = "Time", yscale = log10)
    lines!(time_steps, abundance_true; color = :black, linewidth)
    lines!(time_steps, abundance_lin; color = :black, linewidth, linestyle = :dash)
    isdir("figures") || mkdir("figures")
    save(
        # "/tmp/plot.png",
        "figures/03_data.png",
        fig;
        size = (500, 600),
        px_per_unit = 3,
    )
end
