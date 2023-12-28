using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra
using RareLotkaVolterra
using Statistics

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 30 # Community richness.
n = 2 # Number of communities.
sigma_values = LinRange(0.01, 0.21, 5) # Typical interspecific interaction strength.
n_sigma = length(sigma_values)
z_intensities = LinRange(0.1, 0.9, 5)
initialize_df() = DataFrame(; sigma = [], norm_z = [], error = [])
df_vec = Any[nothing for _ in 1:n_sigma]
Threads.@threads for k in 1:n_sigma
    sigma = sigma_values[k]
    df_temp = initialize_df()
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    coms = create_communities(S, n; create_interaction_matrix)[S]
    for com in coms
        @info "sigma = $sigma"
        Neq = equilibrium_abundance(com)
        for intensity in z_intensities, k in 1:200
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
            end
            push!(df_temp, [sigma, intensity, mean(species_error)])
        end
    end
    df_vec[k] = df_temp
end
df = reduce(vcat, df_vec)
dfp = combine(groupby(df, [:sigma, :norm_z]), :error => mean => :error_mean)
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
