using CairoMakie
using DataFrames
using LinearAlgebra
using RareLotkaVolterra
using Statistics

# Investigate the correlation between the short-term return rates considering or not
# the non-linear terms.
heaviside(x) = x > 0 ? x : 0
"""
Non-linear degree of a perturbation
"""
nonlinear_degree(x, Neq) = sum(heaviside.(abs.(x) - Neq)) / sum(Neq)
function average_return_rate(trajectory, t; observable::Function = norm)
    (log10(observable(trajectory(0))^2) - log10(observable(trajectory(t))^2)) / 2t
end
intensity = 0.99
sigma = 0.2
S = 50 # Community richness.
n = 50 # Number of communities.
n_perturbation = 10_000
t_avg_values = LinRange(10, 500, 5)
tspan = maximum(t_avg_values)
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
coms = create_communities(S, n; create_interaction_matrix)[S]
df = DataFrame(;
    community_id = [],
    prop_x0_rare = [],
    r_nl = [],
    r_lin = [],
    observable = [],
    t_avg = [],
)
for (c_id, com) in enumerate(coms)
    Neq = equilibrium_abundance(com)
    factor = expected_proportional_intensity(Neq)
    Nmin = minimum(Neq)
    observables = Dict(
        :dist => norm, # Distance to equilibrium.
        :mean => mean, # Mean disturbance.
        :rare => (rarest_species(x) = x[argmin(Neq)]), # Distance of the rarest species.
        :rich => (richer_species(x) = x[argmax(Neq)]), # Distance of the most abundant species.
    )
    for i in 1:n_perturbation
        x0 = proportional_perturbation(Neq, factor * intensity, factor)
        # x0 = isotrope_perturbation(Neq, intensity; no_extinction = true)
        trajectory = response(com, x0; tspan)
        for t_avg in t_avg_values
            for (obs_name, obs_fun) in observables
                r_lin = average_return_rate(trajectory.linear, t_avg; observable = obs_fun)
                r_nl =
                    average_return_rate(trajectory.nonlinear, t_avg; observable = obs_fun)
                push!(df, [c_id, x0[argmin(Neq)] / Nmin, r_nl, r_lin, obs_name, t_avg])
            end
        end
    end
    @info "Community $c_id done."
end
transform!(df, [:r_nl, :r_lin] => ByRow((x, y) -> abs(y - x) / abs(x)) => :error_r)
se(x) = std(x) / sqrt(length(x))
ci95(x) = 1.96 * se(x)
pdf = combine(groupby(df, [:observable, :t_avg]), :error_r => mean, :error_r => ci95)

fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel = L"t^\mathrm{avg}",
    ylabel = "Return rate relative error (log)",
    yscale = log10,
)
ln = []
colors = [:black, :red, :blue, :green]
observable = [:dist, :mean, :rare, :rich]
for (obs, color) in zip(observable, colors)
    spdf = pdf[pdf.observable.==obs, :]
    temp = scatterlines!(
        Float64.(spdf.t_avg),
        Float64.(spdf.error_r_mean);
        makercolor = color,
        color = color,
    )
    errorbars!(
        Float64.(spdf.t_avg),
        Float64.(spdf.error_r_mean),
        Float64.(spdf.error_r_ci95);
        whiskerwidth = 3,
        color,
    )
    push!(ln, temp)
end
Legend(fig[0, 1], ln, String.(observable); orientation = :horizontal)
save(
    "figures/error-return-rate-x-proportional.png",
    fig;
    resolution = (450, 350),
    px_per_unit = 3,
)

fig = Figure()
for (i, obs) in enumerate([:dist, :mean, :rare, :rich])
    for (j, t_avg) in enumerate(t_avg_values)
        sdf = df[df.observable.==obs.&&df.t_avg.==t_avg, :]
        x = Float64.(sdf.r_nl)
        y = Float64.(sdf.r_lin)
        rho = round(cor(x, y); digits = 2)
        ax = Axis(
            fig[j, i];
            xlabel = L"R^\mathrm{avg}_{100}",
            ylabel = L"\hat{R}^\mathrm{avg}_{100}",
            title = "Observable $obs, t = $t_avg, cor = $rho",
        )
        strokewidth = 1
        max_deg = maximum(sdf.prop_x0_rare)
        scatter!(
            x,
            y;
            color = Float64.(sdf.prop_x0_rare),
            colorrange = (-1, 1),
            colormap = :delta,
            strokewidth,
        )
        ablines!(0, 1; color = :black)
    end
end
# Colorbar(fig[1, 2]; limits = (-5, 5), colormap = :delta, label = L"x_0^\mathrm{rare}")
save(
    "figures/correlation-return-rate-x-isotrope-observables.png",
    fig;
    resolution = (4 * 450, length(t_avg_values) * 350),
    px_per_unit = 3,
)
