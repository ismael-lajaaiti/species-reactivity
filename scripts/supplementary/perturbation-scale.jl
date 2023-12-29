using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 1 # Number of communities.
sigma = 0.26
intensity_values = rand(Uniform(0, 0.4), 200) # Perturbation intensities.
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
coms = create_communities(S, n; create_interaction_matrix)[S]
df_error = DataFrame(;
    x_norm = Float64[],
    z_norm = Float64[],
    error = Float64[],
    reactivity = Float64[],
)
for c in coms
    Neq = equilibrium_abundance(c)
    A = c.A
    for intensity in intensity_values
        species_error = zeros(S)
        reactivity = zeros(S)
        no_extinction = true
        x0 = isotrope_perturbation(Neq, intensity; no_extinction)
        z0 = x0 ./ Neq
        r = response(c, x0)
        for i in 1:S
            species_error[i] = trajectory_error(
                t -> r.nonlinear(t; idxs = i),
                t -> r.linear(t; idxs = i);
                tspan = (0, 1_000),
            )
            reactivity[i] = sum(sign(z0[i]) * (A[i, :] .* Neq .* z0)) / norm(z0)
        end
        push!(df_error, [norm(x0), norm(z0), mean(species_error), mean(reactivity)])
    end
end

with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Perturbation nonlinearity, $<z_{0, i}^2>_i$",
        ylabel = L"Community linearization error, $<e^{(i)}>_i$",
        xscale = log10,
    )
    ax2 = Axis(
        fig[1, 2];
        xlabel = L"Perturbation intensity, $<x_{0, i}^2>_i$",
        ylabel = L"Perturbation nonlinearity, $<z_{0, i}^2>_i$",
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
    # markersize = 7
    color = :black
    strokecolor = color
    strokewidth = 1
    alpha = 0.6
    scatter!(ax1, (1 / sqrt(S)) .* df_error.z_norm, df_error.error; color)
    limits = extrema(df_error.reactivity)
    scatter!(ax2, (1 / sqrt(S)) .* df_error.x_norm, (1 / sqrt(S)) .* df_error.z_norm; color)
    save("/tmp/plot.png", fig; resolution = (1.1 * 600, 1.1 * 330), px_per_unit = 3)
end
