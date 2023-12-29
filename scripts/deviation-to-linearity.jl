using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics
include("makie-theme.jl")

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 1 # Number of communities.
sigma = 0.26
intensity_values = sqrt(S) .* [0.5, 0.75, 1.0, 1.25, 1.5] # Perturbation intensities.
n_perturbation = 1_000
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
c = create_communities(S, n; create_interaction_matrix)[S][1]
Neq = equilibrium_abundance(c)
A = c.A # - Diagonal(c.A)
remove_diag(A) = A - Diagonal(A)
reactivity(A, Neq, i) = sum((remove_diag(A)[i, :] .* Neq) .^ 2)
reactivity(A, Neq, i, z) = sign(z[i]) * sum(A[i, :] .* Neq .* z) / norm(z)
species_reactivity = [reactivity(A, Neq, i) for i in 1:S]
create_df() = DataFrame(;
    species_id = Float64[],
    abundance = Float64[],
    perturbation_intensity = Float64[],
    perturbation_nonlinearity = Float64[],
    deviation = Float64[],
    reactivity = Float64[],
    reactivity_worst = Float64[],
)
x_set = [
    proportional_perturbation(Neq, maximum(intensity_values), 1, true) for
    _ in 1:n_perturbation
]
df_vector = Any[nothing for intensity_idx in intensity_values]
Threads.@threads for intensity_idx in 1:length(intensity_values)
    intensity = intensity_values[intensity_idx]
    df = create_df()
    for k in 1:n_perturbation
        @info "perturbation $k"
        x0 = (intensity / maximum(intensity_values)) * x_set[k]
        z0 = x0 ./ Neq
        p_intensity = (1 / sqrt(S)) * norm(x0)
        p_nonlinearity = (1 / sqrt(S)) * intensity # (1 / sqrt(S)) * norm(z0)
        r = response(c, x0)
        for i in 1:S
            deviation = trajectory_error(
                t -> r.nonlinear(t; idxs = i),
                t -> r.linear(t; idxs = i);
                tspan = (0, 1_000),
            )
            r0_mean = reactivity(A, Neq, i, z0)
            r0_worst = species_reactivity[i]
            push!(
                df,
                [i, Neq[i], p_intensity, p_nonlinearity, deviation, r0_mean, r0_worst],
            )
        end
    end
    df_vector[intensity_idx] = df
end
df = reduce(vcat, df_vector)
p_df = combine(
    groupby(df, [:species_id, :abundance, :reactivity_worst, :perturbation_nonlinearity]),
    :deviation => mean => :deviation_mean,
    :reactivity => mean => :reactivity_mean,
)

n_x = 500
rand_intensity() = rand(Uniform(0, 1))
x_independent_set = [isotrope_perturbation(Neq, rand_intensity()) for _ in 1:n_x]
z_set = [x ./ Neq for x in x_independent_set]
x_intensity = (1 / sqrt(S)) * norm.(x_independent_set)
z_intensity = (1 / sqrt(S)) * norm.(z_set)

with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2:3] = GridLayout()
    ax1 = Axis(
        fig[1, 2];
        xlabel = L"Perturbation nonlinearity, $<z_{0, i}^2>_i$",
        ylabel = L"Deviation from linearity, $\Delta^{(i)}$",
        # yscale = log10,
    )
    xlim = extrema(p_df.perturbation_nonlinearity)
    width = 0.05 * (xlim[2] - xlim[1])
    jitter() = width .* rand(nrow(p_df))
    scatter!(
        p_df.perturbation_nonlinearity .+ jitter(),
        p_df.deviation_mean;
        color = p_df.reactivity_mean,
    )
    limits = extrema(p_df.reactivity_mean)
    Colorbar(fig[1, 3]; label = L"Mean reactivity, $<R_0^{(i)}>_{\mathbf{z}_0}$", limits)
    ax2 = Axis(
        fig[1, 1];
        xlabel = L"Perturbation intensity, $<x_{0, i}^2>_i$",
        ylabel = L"Perturbation nonlinearity, $<z_{0, i}^2>_i$",
    )
    scatter!(x_intensity, z_intensity; color = :black, alpha = 0.6)
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
        "figures/deviation-to-linearity.png",
        fig;
        resolution = (1.1 * 600, 1.1 * 330),
        px_per_unit = 3,
    )
end

