using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics

S = 30 # Community richness.
n = 1 # Number of communities.
n_perturbation = 1_000
sigma_values = LinRange(0.01, 0.24, 10)
tspan = (0, 10_000)
function initialize_df()
    DataFrame(;
        time_outside_linearity = Float64[],
        maximum_distance = Float64[],
        error = Float64[],
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
        @info "sigma = $sigma, n_min = $(minimum(Neq_relative))"
        intensity = 5
        a = expected_proportional_intensity(Neq)
        for k in 1:n_perturbation
            @info "Perturbation $k"
            x0 = proportional_perturbation(Neq, intensity, 1, true)
            r = response(c, x0; tspan)
            x_true(t) = r.nonlinear(t)
            x_hat(t) = r.linear(t)
            error = fill(0.0, S)
            time = fill(0.0, S)
            tau = fill(0.0, S)
            max_dist = fill(0.0, S)
            for i in 1:S
                tau[i] = transitional_regime_duration(t -> r.nonlinear(t; idxs = i))
                max_dist[i], time[i] =
                    nonlinear_excursion_data(t -> x_true(t)[i], Neq[i]; tspan = (0, tau[i]), delta_t = 0.5)
                error[i] = trajectory_error(
                    t -> r.nonlinear(t; idxs = i),
                    t -> r.linear(t; idxs = i);
                    tspan = (0, tau[i]),
                )
                # time = quadgk(t -> t * abs(r.nonlinear(t; idxs = i)) / Neq[i], 0, 1_000)[1] / 1_000
                # max_dist[i] = quadgk(t -> (abs(r.nonlinear(t; idxs = i)) / Neq[i])^2, 0, time[i])[1] / time[i]
            end
            if all(tau .< 1e5) # Simulation converged.
                for i in 1:S
                    push!(df_tmp, [time[i], max_dist[i], error[i], sigma, Neq_relative[i], i])
                end
            end
        end
    end
    push!(df_vec, df_tmp)
    @info "Sigma = $sigma done."
end
df = reduce(vcat, df_vec)
df_nl_only = df[df.maximum_distance.>0, :]
df_processed = combine(
    groupby(df, [:interaction_strength, :relative_abundance, :species_id]),
    :maximum_distance => mean,
    :error => median => :error_mean,
    :time_outside_linearity => mean,
)
df_nl_only_processed = combine(
    groupby(df_nl_only, [:interaction_strength, :relative_abundance, :species_id]),
    :time_outside_linearity => median => :time_outside_linearity_mean,
    :maximum_distance => median => :maximum_distance_mean,
)

with_theme(publication_theme) do
    fig = Figure()
    # a = fig[1, 1] = GridLayout()
    # b = fig[2, 1] = GridLayout()
    # c = fig[2, 2] = GridLayout()
    color = log10.(df_processed.relative_abundance)
    color_nl = log10.(df_nl_only_processed.relative_abundance)
    colorrange = extrema(color)
    jitter = rand(Uniform(-0.001, 0.001), nrow(df_processed))
    jitter_nl = rand(Uniform(-0.001, 0.001), nrow(df_nl_only_processed))
    ax1 = Axis(
        fig[1, 1:2];
        xlabel = L"Interaction strength, $\sigma$",
        ylabel = L"Linearization error, $e_i$",
        # yscale = log10,
    )
    scatter!(
        df_processed.interaction_strength .+ jitter,
        df_processed.error_mean;
        color,
    )
    Colorbar(fig[1:2, 3]; label = L"log relative abundance, $\log N_i^*$", colorrange)
    ax2 = Axis(
        fig[2, 1];
        xlabel = L"Interaction strength, $\sigma$",
        ylabel = L"Reactive time, $\tau_{i}$",
        yscale = log10,
    )
    scatter!(
        df_nl_only_processed.interaction_strength .+ jitter_nl,
        df_nl_only_processed.time_outside_linearity_mean;
        color = color_nl,
    )
    ax3 = Axis(
        fig[2, 2];
        xlabel = L"Interaction strength, $\sigma$",
        ylabel = L"Max. non-linearity level, $\delta z_\mathrm{max}$",
        # yscale = log10,
        # tellwidth = true,
    )
    scatter!(
        df_nl_only_processed.interaction_strength .+ jitter_nl,
        df_nl_only_processed.maximum_distance_mean;
        color = color_nl,
    )
    save("/tmp/plot.png", fig; resolution = (600, 2 * 320), px_per_unit = 3)
end
