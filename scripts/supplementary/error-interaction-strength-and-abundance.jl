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
n_perturbation = 200
sigma_values = LinRange(0.05, 0.2, 5)
tspan = (0, 100_000)
threshold = 1e-2
function initialize_df()
    DataFrame(;
        deviation_error = Float64[],
        tau = Float64[],
        tau_hat = Float64[],
        interaction_strength = Float64[],
        relative_abundance = Float64[],
        species_id = Int64[],
    )
end
df_vec = DataFrame[]
Threads.@threads for sigma in sigma_values
    df_tmp = initialize_df()
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    coms = create_communities(S, n; create_interaction_matrix)[S]
    for c in coms
        Neq = equilibrium_abundance(c)
        Neq_relative = Neq / sum(Neq)
        intensity = 0.5
        a = expected_proportional_intensity(Neq)
        for k in 1:n_perturbation
            println(k)
            x0 = proportional_perturbation(Neq, intensity, a, true)
            r = response(c, x0; tspan)
            x_true(t) = r.nonlinear(t)
            x_hat(t) = r.linear(t)
            for i in 1:S
                tau = transitional_regime_duration(t -> x_true(t)[i]; threshold)
                tau_hat = transitional_regime_duration(t -> x_hat(t)[i]; threshold)
                e_x = trajectory_error(
                    t -> x_true(t)[i],
                    t -> x_hat(t)[i];
                    tspan = (0, tau),
                )
                push!(df_tmp, [e_x, tau, tau_hat, sigma, Neq_relative[i], i])
            end
        end
    end
    push!(df_vec, df_tmp)
    @info "Sigma = $sigma done."
end
df = reduce(vcat, df_vec)
transform!(df, [:tau, :tau_hat] => ByRow((x, y) -> abs(x - y)) => :e_tau_abs)
transform!(df, [:tau, :tau_hat] => ByRow((x, y) -> abs(x - y) / x) => :e_tau_rel)
df_processed = combine(
    groupby(df, [:interaction_strength, :relative_abundance, :species_id]),
    :deviation_error => mean,
    :e_tau_abs => mean,
    :e_tau_rel => mean,
    :tau => mean,
)

with_theme(publication_theme) do
    fig = Figure()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Interaction strength, $\sigma$",
        ylabel = L"Deviation error, $e_{x, i}$",
    )
    markersize = 8
    color = log10.(df_processed.relative_abundance)
    colorrange = extrema(color)
    jitter = rand(Uniform(-0.001, 0.001), nrow(df_processed))
    scatter!(df_processed.interaction_strength .+ jitter, df_processed.deviation_error_mean; color, markersize)
    Colorbar(fig[1, 2]; label = "log relative abundance", colorrange)
    ax2 = Axis(
        fig[2, 1];
        xlabel = L"Interaction strength, $\sigma$",
        ylabel = L"Return time, $\tau_i$",
        yscale = log10,
    )
    scatter!(
        df_processed.interaction_strength .+ jitter,
        abs.(df_processed.tau_mean);
        color,
        markersize,
    )
    Colorbar(fig[2, 2]; label = "log relative abundance", colorrange)
    ax3 = Axis(
        fig[3, 1];
        xlabel = L"Interaction strength, $\sigma$",
        ylabel = L"Return time abs. error, $e_{\tau, i}$",
    )
    scatter!(
        df_processed.interaction_strength .+ jitter,
        abs.(df_processed.e_tau_abs_mean);
        color,
        markersize,
    )
    ax4 = Axis(
        fig[4, 1];
        xlabel = L"Interaction strength, $\sigma$",
        ylabel = L"Return time rel. error, $e_{\tau, i}$",
    )
    scatter!(
        df_processed.interaction_strength .+ jitter,
        abs.(df_processed.e_tau_rel_mean);
        color,
        markersize,
    )
    save("/tmp/plot.png", fig; resolution = (600, 3 * 320), px_per_unit = 3)
end
