using CairoMakie
using DataFrames
using LinearAlgebra
using RareLotkaVolterra
using Statistics

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n = 50 # Number of communities.
sigma = 0.2 # Typical interspecific interaction strength.
n_perturbation = 100
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
coms = create_communities(S, n; create_interaction_matrix)[S]
intensity_values = range(0.01, 0.1; length = 10)
df = DataFrame(; community_id = [], intensity = [], z_max = [], type = [])
for (c_id, com) in enumerate(coms)
    Neq = equilibrium_abundance(com)
    a = expected_proportional_intensity(Neq)
    for intensity in intensity_values, i in 1:n_perturbation
        x0_iso = isotrope_perturbation(Neq, intensity * sum(Neq))
        x0_prop = proportional_perturbation(Neq, intensity * sum(Neq), a)
        for (type, x0) in Dict(:isotrope => x0_iso, :proportional => x0_prop)
            z0 = x0 ./ Neq
            push!(df, [c_id, intensity, maximum(z0), type])
        end
    end
end
p_df = combine(
    groupby(df, [:intensity, :type]),
    :z_max .=> mean .=> :z_max_mean,
)
transform!(p_df, :intensity => x -> 100x; renamecols = false)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"Intensity ($% N_\mathrm{tot})", ylabel = "Probability")
strokewidth = 1
nl = scatterlines!(Float64.(p_df.intensity), p_df.p_nonlinear)
ext = scatterlines!(Float64.(p_df.intensity), p_df.p_extinction)
Legend(fig[0, :], [nl, ext], ["non-linear", "extinction"]; orientation = :horizontal)
save(
    "figures/probability-x-isotrope-nonlinear.png",
    fig;
    resolution = (450, 350),
    px_per_unit = 3,
)

# Compute the probability of generating a proportional perturbation whose resulting
# trajectory goes outside the linearity domain (even if the proportional perturbation
# is initially within the domain).
"""
Does a given trajectory (of type function) stays within the linear domain?
"""
function is_linear(trajectory, tspan, Neq)
    out, species_order = true, Dict()
    Neq_sorted = sort(Neq)
    prev_species_out = []
    species_out = []
    for t in tspan[1]:1:tspan[2]
        prev_species_out = species_out
        species_out = abs.(trajectory(t) .> Neq)
        if species_out != prev_species_out && any(species_out)
            out = false
            Neq_out = Neq[species_out]
            species_order[t] = [findfirst(==(N), Neq_sorted) for N in Neq_out]
        end
    end
    (out, species_order)
end
t_end = 10_000
tspan = (0, t_end)
n_perturbation = 100
intensity_values = [0.99]
df = DataFrame(; community_id = [], intensity = [], is_linear = [], order = [])
for (c_id, com) in enumerate(coms)
    Neq = equilibrium_abundance(com)
    factor = expected_proportional_intensity(Neq)
    for intensity in intensity_values, i in 1:n_perturbation
        x0 = proportional_perturbation(Neq, factor * intensity, factor)
        r = response(com, x0; tspan)
        islin, order = is_linear(r.nonlinear, tspan, Neq)
        push!(df, [c_id, intensity, islin, order])
    end
    @info "Community $c_id done."
end
p_df = combine(groupby(df, :intensity), :is_linear => mean => :p_linear)

# Probability of escaping trajectory vs. interaction strength 'sigma'
# and reactivity 'r0' (point colors).
"""
Largest initial slope of the community after a pulse perturbation.
"""
function reactivity(com)
    j = jacobian(com)
    0.5 * eigvals(j + j')[end]
end
sigma_values = range(0.17, 0.21; length = 6)
S = 50
n = 100
n_perturbation = 1_000
t_end = 1_000
tspan = (0, t_end)
intensity = 0.99
df = DataFrame(; community_id = [], sigma = [], r0 = [], first_rank = [])
for sigma in sigma_values
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    coms = create_communities(S, n; create_interaction_matrix)[S]
    for (c_id, com) in enumerate(coms)
        Neq = equilibrium_abundance(com)
        factor = expected_proportional_intensity(Neq)
        r0 = -reactivity(com)
        for _ in 1:n_perturbation
            x0 = proportional_perturbation(Neq, factor * intensity, factor)
            r = response(com, x0; tspan)
            islin, order = is_linear(r.nonlinear, tspan, Neq)
            first_rank = isempty(order) ? nothing : mean(order[minimum(keys(order))])
            push!(df, [c_id, sigma, r0, first_rank])
        end
    end
    @info "Sigma $sigma done."
end
pdf = combine(
    groupby(df, :sigma),
    :first_rank => (x -> replace(x, nothing => 0) |> mean) => :p_escape,
    :r0 => mean,
)
spdf = pdf[pdf.p_escape.!=0, :]

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"\sigma", ylabel = L"{p}_\mathrm{escape}")
strokewidth = 1
scatter!(
    Float64.(spdf.sigma),
    log10.(Float64.(spdf.p_escape));
    color = Float64.(spdf.r0_mean),
    colorrange = extrema(spdf.r0_mean),
    colormap = :viridis,
    strokewidth,
)
Colorbar(
    fig[1, 2];
    limits = extrema(spdf.r0_mean),
    colormap = :viridis,
    label = L"R^\mathrm{inst}_{0}",
)
ax = Axis(fig[1, 3]; xlabel = "Abundance rank", ylabel = "Frequency", xticks = [1, 2, 3])
hist!(
    Int64.(df[.!isnothing.(df.first_rank), :].first_rank);
    bins = 0.5 .+ [0, 1, 2, 3],
    color = :gray,
    strokecolor = :black,
    normalization = :pdf,
)
save(
    "figures/probability-escaping-trajectories.png",
    fig;
    resolution = (700, 350),
    px_per_unit = 3,
)
