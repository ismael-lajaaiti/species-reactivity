using DataFrames
using Distributions
using LinearAlgebra
using QuadGK
using RareLotkaVolterra
using Statistics
include("makie-theme.jl")

function find_max(x; tspan = (0, 1_000))
    t0, tend = tspan
    time_steps = t0:1:tend
    maximum(x.(time_steps))
end

median_reactivity(jacobian) = mean(real.(eigvals(jacobian)))

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 1 # Number of communities.
n_perturbation = 1_000
perturbation_intensity = 5
sigma_values = collect(0.05:0.05:0.25)
df = DataFrame(;
    sigma = Float64[],
    reactivity = Float64[],
    mean_nonlinear = Float64[],
    max_nonlinear = Float64[],
)
for sigma in sigma_values
    @info "sigma = $sigma"
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
    Neq = equilibrium_abundance(com)
    rare = argmin(Neq)
    pseudo_jacobian = com.A * Diagonal(Neq)
    reac = median_reactivity(pseudo_jacobian)
    for k in 1:n_perturbation
        @info "Perturbation $k"
        x0 = proportional_perturbation(Neq, perturbation_intensity, 1, true)
        r = response(com, x0)
        df_tmp = DataFrame(; mean_nl = Float64[], max_nl = Float64[], tau = Float64[])
        for i in 1:S
            tau = transitional_regime_duration(t -> r.nonlinear(t; idxs = i))
            nl(t) = abs(r.nonlinear(t; idxs = i)) / Neq[i]
            mean_nl = quadgk(nl, 0, tau)[1] / tau
            max_nl = find_max(nl; tspan = (0, tau))
            push!(df_tmp, [mean_nl, max_nl, tau])
        end
        if maximum(df_tmp.tau) < 10_000
            mean_nl = df_tmp.mean_nl[rare]
            max_nl = df_tmp.max_nl[rare]
            push!(df, [sigma, reac, mean_nl, max_nl])
        end
    end
end
pdf = combine(
    groupby(df, [:sigma, :reactivity]),
    :mean_nonlinear => mean => :mean_nonlinear_mean,
    :max_nonlinear => mean => :max_nonlinear_mean,
)

with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Reactivity, $r(A)$",
        ylabel = L"Community mean non-linear level, $\mathcal{Z}_\mathrm{mean}$",
        yscale = log10,
    )
    scatter!(pdf.reactivity, pdf.mean_nonlinear_mean; color = pdf.sigma)
    ax2 = Axis(
        fig[1, 2];
        xlabel = L"Reactivity, $r(A)$",
        ylabel = L"Maximum non-linear level, $\max_t |z_i|$",
    )
    scatter!(pdf.reactivity, pdf.max_nonlinear_mean; color = pdf.sigma)
    limits = extrema(pdf.sigma)
    Colorbar(fig[1, 3]; label = L"Interaction strength, $\sigma", limits)
    for (layout, label) in zip([a, b], ["A", "B"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    save("/tmp/plot.png", fig; resolution = (620, 320), px_per_unit = 3)
end


