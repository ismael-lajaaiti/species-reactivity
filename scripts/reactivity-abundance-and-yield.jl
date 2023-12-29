include("makie-theme.jl")
using DataFrames
using Distributions
using GLM
using LinearAlgebra
using QuadGK
using RareLotkaVolterra
using Statistics

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 1 # Number of communities.
sigma = 0.26
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
A = com.A
yield = equilibrium_abundance(com)
remove_diagonal(A) = A - Diagonal(A)
get_reactivity(A, Neq, i) = sqrt(sum((remove_diagonal(A)[i, :] .* Neq) .^ 2))
reactivity = [get_reactivity(A, yield, i) for i in 1:S]
K = rand(LogNormal(0, 0.5), S) ./ yield
abundance = yield .* K
df = DataFrame(; reactivity, abundance, yield)
ols_yield = lm(@formula(reactivity ~ yield), df)
ols_abundance = lm(@formula(reactivity ~ abundance), df)

with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    for (layout, label) in zip([a, b], ["A", "B"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Abundance, $N_i^*$",
        ylabel = L"Reactivity, $R_0^{(i)}$",
    )
    color = :black
    alpha = 1
    scatter!(abundance, reactivity; color, alpha)
    # Linear regression.
    a, b = coef(ols_abundance)
    r2_abundance = round(r2(ols_abundance); digits = 2)
    xlims = [minimum(abundance), maximum(abundance)]
    ylims = a .+ b .* xlims
    lines!(xlims, ylims; color = :grey, label = L"r^2 = %$r2_abundance")
    axislegend()
    ax2 = Axis(
        fig[1, 2];
        xlabel = L"Relative yield, $\eta_i^*$",
    )
    scatter!(yield, reactivity; color, alpha)
    # Linear regression.
    a, b = coef(ols_yield)
    r2_yield = round(r2(ols_yield); digits = 2)
    xlims = [minimum(yield), maximum(yield)]
    ylims = a .+ b .* xlims
    lines!(xlims, ylims; color = :grey, label = L"r^2 = %$r2_yield")
    axislegend()
    save("figures/reactivity-abundance-and-yield.png", fig; resolution = (620, 320), px_per_unit = 3)
end

