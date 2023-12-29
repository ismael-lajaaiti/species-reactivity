using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using SpeciesReactivity
using Statistics


# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n = 1 # Number of communities.
n_perturbation = 50
sigma = 0.1 # Typical interspecific interaction strength.
time_steps = LinRange(10, 500, 7)
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
coms = create_communities(S, n; create_interaction_matrix)[S]
df = DataFrame(; error = Float64[], nb_species = Float64[], x_sign = Symbol[])
for c in coms
    Neq = equilibrium_abundance(c)
    intensity = 0.1
    a = expected_proportional_intensity(Neq)
    for k in 1:n_perturbation
        @info "Perturbation: $k"
        x0 = abs.(proportional_perturbation(Neq, intensity, a, true))
        for x_sign in [:positive, :negative]
            x_sign == :negative && (x0 = -x0)
            r = response(c, x0)
            for n_species in 1:S
                sampled_species = sample(1:S, n_species; replace = false)
                x_true(t) = sum(r.nonlinear(t)[sampled_species])
                x_hat(t) = sum(r.linear(t)[sampled_species])
                e = trajectory_error(x_true, x_hat; tspan = (0, 500))
                push!(df, [e, n_species, x_sign])
            end
        end
    end
end
pdf = combine(groupby(df, [:nb_species, :x_sign]), :error => mean)

# Example of estimated vs. true recovery trajectory of species abundances.
# For most abundant and least abundant species.
time_steps_ex = LinRange(0, 500, 1_000)
df_example = DataFrame(;
    t = Float64[],
    x_true = Float64[],
    x_hat = Float64[],
    collectivity = Symbol[],
)
c = coms[1]
Neq = equilibrium_abundance(c)
abundance_rank = invperm(sortperm(Neq)) # From least abundant (1) to most abundant (S).
single_sp = abundance_rank[div(S, 2)]
a = expected_proportional_intensity(Neq)
x0 = abs.(proportional_perturbation(Neq, 1, a, true))
r = response(c, x0)
e_com = trajectory_error(t -> norm(r.nonlinear(t)), t -> norm(r.linear(t)); tspan = (0, 100))
e_single = trajectory_error(t -> norm(r.nonlinear(t)[single_sp]), t -> norm(r.linear(t)[single_sp]); tspan = (0, 100))
for t in time_steps_ex
    for traj in [:single_sp, :community]
        tau = traj == :single_sp ? t : t / 50
        x_true = r.nonlinear(tau)
        x_hat = r.linear(tau)
        if traj == :community
            traj_true = sum(x_true)
            traj_hat = sum(x_hat)
        elseif traj == :single_sp
            traj_true = sum(x_true[single_sp])
            traj_hat = sum(x_hat[single_sp])
        end
        push!(df_example, [tau, traj_true, traj_hat, traj])
    end
end

with_theme(publication_theme) do
    fig = Figure()
    a = fig[1, 1:2] = GridLayout()
    row2 = fig[2, 1:2] = GridLayout()
    b = row2[1, 1]
    c = row2[1, 2]
    ax1 = Axis(
        a[1, 1];
        xlabel = L"Size of community subset, $s$",
        ylabel = L"Error per community subset, $e_{<\mathbf{x}>}$",
    )
    markersize = 10
    pdf_pos = pdf[pdf.x_sign .== :positive, :]
    pdf_neg = pdf[pdf.x_sign .== :negative, :]
    scatter!(pdf_pos.nb_species, pdf_pos.error_mean; color = :blue, markersize, label = L"x_i(0) \geq 0")
    scatter!(pdf_neg.nb_species, pdf_neg.error_mean; color = :red, markersize, label = L"x_i(0) \leq 0")
    axislegend()
    ax2 = Axis(
        b;
        xlabel = L"Time, $t$",
        ylabel = L"Averaged deviation, $<\mathbf{x}>$",
        title = L"\text{Single species}",
    )
    df_example_min = df_example[df_example.collectivity.==:single_sp, :]
    color = :black
    lines!(df_example_min.t, df_example_min.x_true; color)
    lines!(df_example_min.t, df_example_min.x_hat; color, linestyle = :dash)
    ax3 = Axis(c; xlabel = L"Time, $t$", title = L"\text{Community}")
    df_example_max = df_example[df_example.collectivity.==:community, :]
    lines!(df_example_max.t, df_example_max.x_true; color, label = L"<\mathbf{x}>")
    lines!(
        df_example_max.t,
        df_example_max.x_hat;
        color,
        linestyle = :dash,
        label = L"<\mathbf{\hat{x}}>",
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

