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
n_perturbation = 1_000
perturbation_intensity = 5
sigma_values = collect(0.05:0.05:0.25)
df = DataFrame(;
    abundance = Float64[],
    return_time = Float64[],
    mean_nonlinear = Float64[],
    z0_i_squared = Float64[],
    sigma = Float64[],
)
for sigma in sigma_values
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
    Neq = equilibrium_abundance(com)
    for k in 1:n_perturbation
        @info "Perturbation $k"
        x0 = proportional_perturbation(Neq, perturbation_intensity, 1, true)
        r = response(com, x0)
        for i in 1:S
            return_time = transitional_regime_duration(t -> r.nonlinear(t; idxs = i))
            tau = 1_000
            nl(t) = abs(r.linear(t; idxs = i)) / Neq[i]
            mean_nl = quadgk(t -> nl(t)^2, 0, tau)[1] # / tau
            # error = trajectory_error(
            #     t -> r.nonlinear(t; idxs = i),
            #     t -> r.linear(t; idxs = i);
            #     tspan = (0, tau),
            # )
            push!(df, [Neq[i], return_time, mean_nl, (x0[i] / Neq[i])^2, sigma])
        end
    end
end
df = df[df.return_time.<10_000, :]
pdf = combine(
    groupby(df, :abundance),
    :error => mean,
    :mean_nonlinear => mean,
    :max_nonlinear => mean,
    :z0_i_squared => mean,
)
transform!(
    pdf,
    [:z0_i_squared_mean, :abundance] => ByRow((x, y) -> x / (2 * y)) => :expected,
)
sort!(pdf, :abundance)

with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Species abundance, $N_i^*$",
        ylabel = L"Mean non-linear level, $<|z_i|>_t$",
        # yscale = log10,
    )
    scatter!(pdf.abundance, pdf.mean_nonlinear_mean; color = pdf.error_mean)
    lines!(pdf.abundance, pdf.expected; color = :black)
    ax2 = Axis(
        fig[1, 2];
        xlabel = L"Species abundance, $N_i^*$",
        ylabel = L"Maximum non-linear level, $\max_t |z_i|$",
    )
    scatter!(pdf.abundance, pdf.max_nonlinear_mean; color = pdf.error_mean)
    limits = extrema(pdf.error_mean)
    Colorbar(fig[1, 3]; label = L"Linearization error, $e_\mathrm{lin}$", limits)
    for (layout, label) in zip([a, b], ["A", "B"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    save("figures/non-linear-level-lin.png", fig; resolution = (620, 320), px_per_unit = 3)
end

