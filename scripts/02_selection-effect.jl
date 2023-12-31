using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using GLM
using LinearAlgebra
using QuadGK
using Random
using SpeciesReactivity
using Statistics

Random.seed!(112)

include("makie-theme.jl") # For reproducibility.

S = 30
n = 1
sigma = 0.25
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
eta = equilibrium_abundance(com)
color_scheme = ColorSchemes.seaborn_pastel

K = eta .* (1 .+ rand(Normal(0, 0.2), S))
N = K .* eta
K = rand(Normal(1, 0.2), S)
K = -eta .* (1 .+ rand(Normal(0, 0.3), S)) .+ 1.5 * maximum(eta)


# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 1 # Number of communities.
sigma = 0.26
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
A = com.A
eta = equilibrium_abundance(com)
reactivity = [get_reactivity(A, eta, i) for i in 1:S]

K1 = (1 .+ rand(Normal(0, 0.1), S)) ./ eta .^ 1.25
# K1 = (10 .+ rand(Normal(0, 0.1), S) .- 10 * eta) ./ eta
K2 = (1 .+ rand(Normal(0, 0.1), S)) ./ eta
# K3 = (1 .+ rand(Normal(0, 0.1), S))
K3 = (1 .+ rand(Normal(0, 0.1), S)) .* eta
K_vec = [K1, K2, K3]
idx = argmin(eta)
x0 = remove_diagonal(A)[idx, :] .* eta
r = SpeciesReactivity.response(com, x0)
time_steps = collect(0:0.1:50)
eta_true = [abs.(r.nonlinear(t)) for t in time_steps]
eta_lin = [abs.(r.linear(t)) for t in time_steps]

with_theme(p_theme) do
    fig = Figure()
    strokewidth = 0.5
    markersize = 9
    title_vec = [L"\mathrm{Cov}(N, η)<0", L"\mathrm{Cov}(N, η)=0", L"\mathrm{Cov}(N, η)>0"]
    color = :grey
    for (i, K) in enumerate(K_vec)
        title = title_vec[i]
        abundance = eta .* K
        relative_abundance = abundance / sum(abundance)
        # Abundance vs. relative yield.
        ylabel = i == 1 ? "Abundance (N)" : ""
        ax = Axis(fig[3, i]; xlabel = "Relative yield (η)", ylabel, title)
        scatter!(ax, eta, relative_abundance; color, strokewidth, markersize)
        # Reactivity vs. abundance: scatter plot.
        ylabel = i == 1 ? "Reactivity" : ""
        ax2 = Axis(fig[4, i]; xlabel = "Abundance", ylabel)
        scatter!(ax2, relative_abundance, reactivity; color, strokewidth, markersize)
        # Reactivity vs. abundance: linear fit.
        df = DataFrame(; reactivity, relative_abundance, eta)
        ols = lm(@formula(reactivity ~ relative_abundance), df)
        p_value = ftest(ols.model).pval
        a, b = coef(ols)
        xlims = [minimum(relative_abundance), maximum(relative_abundance)]
        ylims = a .+ b .* xlims
        linestyle = p_value > 0.05 ? :dot : :solid
        lines!(xlims, ylims; color = :skyblue, linestyle, linewidth = 2)
        # Total abundance vs. time.
        abundance_true = [sum(eta_t .* K) for eta_t in eta_true] / sum(abundance)
        abundance_lin = [sum(eta_t .* K) for eta_t in eta_lin] / sum(abundance)
        # abundance_true = [eta_t[idx] * K[idx] for eta_t in eta_true]
        # abundance_lin = [eta_t[idx] * K[idx] for eta_t in eta_lin]
        ylabel = i == 1 ? "Δ tot. abundance" : ""
        ax = Axis(fig[5, i]; xlabel = "Time", ylabel, yscale = log10)
        lines!(time_steps, abundance_true; color = :black, label = "true", linewidth = 2)
        lines!(
            time_steps,
            abundance_lin;
            color = :black,
            label = "prediction",
            linewidth = 2,
            linestyle = :dash,
        )
        i == 3 && axislegend(ax; labelsize = 10)
    end
    ax = Axis(fig[1, 2]; xlabel = "Relative yield", ylabel = "Reactivity")
    scatter!(eta, reactivity; color, strokewidth, markersize)
    eta_linrange = LinRange(0, maximum(eta), 100)
    exp_reactivity_full = sqrt.(expected_reactivity_squared.(eta_linrange, Ref(com)))
    lines!(
        eta_linrange,
        exp_reactivity_full;
        color = :coral,
        label = "expectation",
        alpha = 0.7,
        linewidth = 2,
    )
    # axislegend(; position = :lb)
    Label(fig[2, 1:2], "Selection effect < 0"; tellwidth = false, fontsize = 20)
    Label(fig[2, 3], "Selection effect ≥ 0"; tellwidth = false, fontsize = 20)
    isdir("figures") || mkdir("figures")
    save(
        "figures/02_selection-effect.png",
        # "/tmp/plot.png",
        fig;
        size = 1.1 .* (620, 620),
        px_per_unit = 3,
    )
end
