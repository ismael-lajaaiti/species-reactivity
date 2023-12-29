using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics

function find_max(x; tspan = (0, 1_000))
    t0, tend = tspan
    time_steps = t0:1:tend
    maximum(x.(time_steps))
end

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 1 # Number of communities.
intensity = 0.1
delta_t = 1e-3
sigma_values = LinRange(0.2, 0.21, 5) # Typical interspecific interaction strength.
df = DataFrame(; N = Float64[], init_slope = Float64[])
for sigma in [0.2]
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    com = create_communities(S, n; create_interaction_matrix)[S][1]
    Neq = equilibrium_abundance(com)
    for k in 1:1_000
        @info "Perturbation $k"
        no_extinction = true
        x0 = proportional_perturbation(Neq, intensity, 1, no_extinction)
        z0 = x0 ./ Neq
        r = response(com, x0)
        for i in 1:S
            z(t) = r.linear(t; idxs = i) / Neq[i]
            initial_slope = (abs(z(delta_t)) - abs(z(0))) / delta_t
            # tau = transitional_regime_duration(t -> r.nonlinear(t; idxs = i))
            # zlin = quadgk(t -> abs(r.linear(t; idxs = i)) / Neq[i], 0, tau; rtol = 1e-2)[1] / tau
            # znl = quadgk(t -> abs(r.linear(t; idxs = i)) / Neq[i], 0, tau; rtol = 1e-2)[1] / tau
            # zlin = find_max(t -> abs(r.linear(t; idxs = i)) / Neq[i]; tspan = (0, tau))
            # maxi = find_max(t -> abs(r.nonlinear(t; idxs = i)) / Neq[i]; tspan = (0, tau))
            push!(df, [Neq[i], initial_slope])
        end
    end
end
pdf = combine(groupby(df, :N), :init_slope => mean)
scatter(pdf.N, pdf.init_slope_mean)

n_sigma = length(sigma_values)
intensity = 1
initialize_df() =
    DataFrame(; sigma = [], reactivity = [], norm_z = [], rarity = [], error = [])
df_vec = Any[nothing for _ in 1:n_sigma]
Threads.@threads for k in 1:n_sigma
    sigma = sigma_values[k]
    df_temp = initialize_df()
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    coms = create_communities(S, n; create_interaction_matrix)[S]
    for com in coms
        @info "sigma = $sigma"
        Neq = equilibrium_abundance(com)
        reactivity = nonlinear_reactivity(com)
        for k in 1:200
            @info "Perturbation $k"
            species_error = zeros(S)
            no_extinction = true
            x0 = proportional_perturbation(Neq, intensity, 1, no_extinction)
            z0 = x0 ./ Neq
            r = response(com, x0)
            for i in 1:S
                species_error[i] = trajectory_error(
                    t -> r.nonlinear(t; idxs = i),
                    t -> r.linear(t; idxs = i);
                    tspan = (0, 1_000),
                )
                rarity = Neq[i] / sum(Neq)
                push!(df_temp, [sigma, reactivity, intensity, rarity, species_error[i]])
            end
        end
    end
    df_vec[k] = df_temp
end
df = reduce(vcat, df_vec)
dfp = combine(groupby(df, [:sigma, :reactivity, :rarity]), :error => mean => :error_mean)

get_row(df, sigma, norm_z) = findfirst(==(1), df.sigma .== sigma .&& df.norm_z .== norm_z)
error_matrix = [
    dfp[get_row(dfp, sigma, norm_z), :error_mean] for sigma in sigma_values,
    norm_z in z_intensities
]

with_theme(publication_theme) do
    fig = Figure()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Interaction strength, $\sigma$",
        ylabel = L"Perturbation non-linearity level, $||\mathbf{z}(0)||$",
    )
    hm = heatmap!(sigma_values, z_intensities, error_matrix)
    Colorbar(fig[1, 2], hm; label = L"Linearization error, $e_\text{tot}")
    save("/tmp/plot.png", fig; resolution = (600, 320), px_per_unit = 3)
end
