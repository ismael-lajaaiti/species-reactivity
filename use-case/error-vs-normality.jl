using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics


# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n = 50 # Number of communities.
n_perturbation = 50
sigma_min = 0.1 # Typical interspecific interaction strength.
sigma_max = 0.21 # Typical interspecific interaction strength.
create_interaction_matrix(S) =
    random_competition_matrix(S, rand(Uniform(sigma_min, sigma_max)))
coms = create_communities(S, n; create_interaction_matrix)[S]
df = DataFrame(;
    c_id = Int64[],
    error = Float64[],
    complexity = Float64[],
    departure_from_normality = Float64[],
)
for (c_id, c) in enumerate(coms)
    @info "Community: $c_id"
    Neq = equilibrium_abundance(c)
    cplx = complexity(c)
    dep_norm = departure_from_normality(c)
    intensity = 0.8
    a = expected_proportional_intensity(Neq)
    for k in 1:n_perturbation
        x0 = abs.(proportional_perturbation(Neq, intensity, a, true))
        r = response(c, x0)
        x_true(t) = r.nonlinear(t)
        x_hat(t) = r.linear(t)
        e = trajectory_error(x_true, x_hat; tspan = (0, 500))
        push!(df, [c_id, e, cplx, dep_norm])
    end
end
pdf = combine(groupby(df, [:c_id, :complexity, :departure_from_normality]), :error => mean)

with_theme(publication_theme) do
    fig = Figure()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Complexity, $||A||_F$",
        ylabel = L"Linearization error, $e$",
    )
    markersize = 10
    scatter!(
        pdf.complexity,
        pdf.error_mean;
        color = :black,
        markersize,
    )
    save("/tmp/plot.png", fig; resolution = (600, 320), px_per_unit = 3)
end

