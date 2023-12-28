using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n = 1 # Number of communities.
n_perturbation = 100
sigma = 0.21 # Typical interspecific interaction strength.
df = DataFrame(; error = Float64[], relative_abundance = Float64[], abundance_rank = Float64[])
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
coms = create_communities(S, n; create_interaction_matrix)[S]
for c in coms
    Neq = equilibrium_abundance(c)
    Neq_relative = Neq / sum(Neq)
    abundance_rank = invperm(sortperm(Neq)) # From least abundant (1) to most abundant (S).
    intensity = 0.8
    a = expected_proportional_intensity(Neq)
    for k in 1:n_perturbation
        println(k)
        x0 = proportional_perturbation(Neq, intensity, a, true)
        r = response(c, x0)
        x_true(t) = r.nonlinear(t)
        x_hat(t) = r.linear(t)
        for i in 1:S
            e = trajectory_error(t -> x_true(t)[i], t -> x_hat(t)[i]; tspan = (0, 500))
            push!(df, [e, Neq_relative[i], abundance_rank[i]])
        end
    end
end
pdf = combine(groupby(df, [:relative_abundance]), :error => mean)

# Example of estimated vs. true recovery trajectory of species abundances.
# For most abundant and least abundant species.
time_steps_ex = LinRange(0, 1_000, 1_000)
df_example =
    DataFrame(; t = Float64[], N_true = Float64[], N_hat = Float64[], species = Symbol[])
c = coms[1]
Neq = equilibrium_abundance(c)
a = expected_proportional_intensity(Neq)
x0 = proportional_perturbation(Neq, 1, a, true)
r = response(c, x0)
index_dict = Dict(:min => argmin(Neq), :max => argmax(Neq))
Nmin = minimum(Neq)
Nmax = maximum(Neq)
# push!(df_example, [-2 / Nmin, Nmin, Nmin, :min])
# push!(df_example, [-2 / Nmax, Nmax, Nmax, :max])
# push!(df_example, [-0.1 / Nmin, Nmin, Nmin, :min])
# push!(df_example, [-0.1 / Nmax, Nmax, Nmax, :max])
e_min = trajectory_error(t -> r.nonlinear(t)[index_dict[:min]], t -> r.linear(t)[index_dict[:min]], tspan = (0, 1_000))
e_max = trajectory_error(t -> r.nonlinear(t)[index_dict[:max]], t -> r.linear(t)[index_dict[:max]], tspan = (0, 1_000))
for t in time_steps_ex
    for species in [:min, :max]
        index = index_dict[species]
        tau = t
        x_true = r.nonlinear(tau)
        x_hat = r.linear(tau)
        push!(df_example, [tau, x_true[index], x_hat[index], species])
    end
end

with_theme(publication_theme) do
    fig = Figure()
    a = fig[1, 1:2] = GridLayout()
    colorbar_grid = fig[1, 2] = GridLayout()
    row2 = fig[2, 1:2] = GridLayout()
    b = row2[1, 1]
    c = row2[1, 2]
    ax1 = Axis(
        a[1, 1];
        xlabel = L"Relative abundance $n_i^*$",
        ylabel = L"Error per species, $e_i$",
    )
    color = :black
    markersize = 10
    scatter!(
        pdf.relative_abundance,
        pdf.error_mean;
        color,
        markersize
    )
    ax2 = Axis(
        b;
        xlabel = L"Time, $t$",
        ylabel = L"Deviation from equilibrium, $x$",
        title = L"\text{Rare}",
    )
    df_example_min = df_example[df_example.species.==:min, :]
    color = :black
    lines!(df_example_min.t, df_example_min.N_true; color)
    lines!(df_example_min.t, df_example_min.N_hat; color, linestyle = :dash)
    ax3 = Axis(c; xlabel = L"Time, $t$", title = L"\text{Abundant}")
    df_example_max = df_example[df_example.species.==:max, :]
    lines!(df_example_max.t, df_example_max.N_true; color, label = L"x")
    lines!(
        df_example_max.t,
        df_example_max.N_hat;
        color,
        linestyle = :dash,
        label = L"\hat{x}",
    )
    axislegend()
    for (layout, label) in zip([a, b, c], ["A", "B", "C"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    save("/tmp/plot.png", fig; resolution = (600, 2 * 320), px_per_unit = 3)
end
