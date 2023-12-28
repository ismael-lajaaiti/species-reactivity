using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics
include("makie-theme.jl")

"""
Average return rate between 0 and t.
"""
return_rate(x, t) = (log(x(0)^2) - log(x(t)^2)) / t

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 1 # Number of communities.
n_perturbation = 5_000
perturbation_intensity = sqrt(S) # <z_i^2> = 1.
sigma = 0.26
tau = 100
df = DataFrame(;
    abundance = Float64[],
    return_rate_true = Float64[],
    return_rate_hat = Float64[],
    deviation = Float64[],
)
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
Neq = equilibrium_abundance(com)
for k in 1:n_perturbation
    @info "Perturbation $k"
    x0 = proportional_perturbation(Neq, perturbation_intensity, 1, true)
    r = response(com, x0)
    z0 = abs.(x0) ./ Neq
    for i in 1:S
        r_true = return_rate(t -> r.nonlinear(t; idxs = i), tau)
        r_hat = return_rate(t -> r.linear(t; idxs = i), tau)
        deviation = trajectory_error(
            t -> r.nonlinear(t; idxs = i),
            t -> r.linear(t; idxs = i);
            tspan = (0, 10_000),
        )
        push!(df, [Neq[i], r_true, r_hat, deviation])
    end
end
df = df[abs.(df.return_rate_true).>1e-3, :]
transform!(
    df,
    [:return_rate_true, :return_rate_hat] => ByRow((x, y) -> abs(x - y)) => :error_abs,
    [:return_rate_true, :return_rate_hat] =>
        ByRow((x, y) -> 100 * abs(x - y) / abs(x)) => :error_rel,
)
pdf = combine(
    groupby(df, :abundance),
    :return_rate_true => mean,
    :return_rate_hat => mean,
    :error_abs => mean => :error_abs_mean,
    :error_rel => median => :error_rel_mean,
    :deviation => mean,
)

with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    colors = [:red, :green, :blue]
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Relative yield, $\eta_i^*$",
        ylabel = L"Error on the return rate, $e_R$ (%)",
        # yscale = log10,
    )
    color = pdf.deviation_mean
    colorrange = extrema(pdf.deviation_mean)
    scatter!(pdf.abundance, pdf.error_rel_mean; color)
    Colorbar(fig[1, 2]; label = L"Deviation from linearity, $\Delta^{(i)}$", colorrange)
    save("figures/error-return-rate.png", fig; resolution = (520, 320), px_per_unit = 3)
end
