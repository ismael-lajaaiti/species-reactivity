using DataFrames
using Distributions
using LinearAlgebra
using QuadGK
using RareLotkaVolterra
using Statistics

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 1 # Number of communities.
n_perturbation = 500
perturbation_intensity = 0.5
sigma = 0.25
delta_t = 0.1
df = DataFrame(;
    species_id = Int64[],
    abundance = Float64[],
    return_time = Float64[],
    initial_derivative = Float64[],
)
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
A = com.A
A_without_diag = A - Diagonal(A)
Neq = equilibrium_abundance(com)
for k in 1:S
    @info "Perturbation $k"
    # x0 = proportional_perturbation(Neq, perturbation_intensity, 1, true)
    x0 = [A_without_diag[k, j] * Neq[j]^2 for j in 1:S]
    r = response(com, x0)
    norm_z0 = norm(x0 ./ Neq)
    for i in 1:S
        return_time = transitional_regime_duration(t -> r.linear(t; idxs = i); threshold = 0.9)
        tau = 1_000
        nl(t) = abs(r.linear(t; idxs = i)) / Neq[i]
        initial_derivative = (nl(delta_t) - nl(0)) / (delta_t * norm_z0)
        push!(df, [i, Neq[i], return_time, initial_derivative])
    end
end
# df = df[df.return_time.<10_000, :]
pdf = combine(
    groupby(df, [:species_id, :abundance]),
    :initial_derivative => mean,
    :initial_derivative => maximum,
)
transform!(
    pdf,
    :species_id => ByRow(i -> sqrt(sum((A_without_diag[i, :] .* Neq) .^ 2))) => :reactivity,
)
sort!(pdf, :abundance)

with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Species reactivity, $R_{0, i}$",
        ylabel = L"Initial derivative, $\dot{z}_i(0)$",
        # yscale = log10,
    )
    scatter!(pdf.reactivity, pdf.initial_derivative_maximum; color = :black)
    xmin, xmax = extrema(pdf.reactivity)
    lines!(LinRange(xmin, xmax, 10), LinRange(xmin, xmax, 10))
    save(
        "figures/reactivity-check.png",
        fig;
        resolution = (620, 320),
        px_per_unit = 3,
    )
end

