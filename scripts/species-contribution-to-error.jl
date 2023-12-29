using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 20 # Number of communities.
n_perturbation = 500
sigma = 0.2 # Typical interspecific interaction strength.
time_steps = [1, 10, 100]
df = DataFrame(; error = Float64[], z = Float64[], n = Float64[], t = Float64[])
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
coms = create_communities(S, n; create_interaction_matrix)[S]
for c in coms
    Neq = equilibrium_abundance(c)
    intensity = 1
    a = expected_proportional_intensity(Neq)
    for _ in 1:n_perturbation
        x0 = proportional_perturbation(Neq, intensity, a, true)
        r = response(c, x0)
        for t in time_steps
            error = abs.(r.linear(t) - r.nonlinear(t))
            error /= norm(error)
            z = r.nonlinear(t) ./ Neq
            for i in 1:S
                push!(df, [error[i] / Neq[i], abs(z[i]), Neq[i] / sum(Neq), t])
            end
        end
    end
end
pdf = combine(groupby(df, [:t, :n]), :error => mean, :z => mean)

colors = [:red, :green, :blue]
with_theme(publication_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Relative abundance, $n^*_i$",
        ylabel = L"Species contribution to error, $\tilde{e}_i$",
        yscale = log10,
    )
    markersize = 7
    strokewidth = 1
    alpha = 0.5
    for (t, color) in zip(time_steps, colors)
        spdf = pdf[pdf.t.==t, :]
        label = L"t = %$t"
        scatter!(spdf.n, spdf.error_mean; markersize, color, label, alpha)
    end
    labelsize = 10
    titlesize = 10
    axislegend()
    ax2 = Axis(
        fig[1, 2];
        xlabel = L"Relative abundance, $n^*_i$",
        ylabel = L"Scaled distance to equilibrium, $|z_i|$",
        yscale = log10,
    )
    markersize = 7
    strokewidth = 1
    alpha = 0.5
    for (t, color) in zip(time_steps, colors)
        spdf = pdf[pdf.t.==t, :]
        label = L"t = %$t"
        scatter!(spdf.n, spdf.z_mean; markersize, color, label, alpha)
    end
    labelsize = 10
    titlesize = 10
    for (layout, label) in zip([a, b], ["A", "B"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    save("/tmp/plot.png", fig; resolution = (600, 320), px_per_unit = 3)
end
