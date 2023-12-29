using CairoMakie
using DataFrames
using LinearAlgebra
using RareLotkaVolterra
using Statistics

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
first_rank_vec = df[.!isnothing.(df.first_rank), :].first_rank
ax = Axis(fig[1, 3]; xlabel = "Abundance rank", ylabel = "Frequency", xticks = [i for i in 1:maximum(first_rank_vec)])
hist!(
    Int64.(first_rank_vec);
    bins = 0.5 .+ [i for i in 0:maximum(first_rank_vec)],
    color = :gray,
    strokewidth = 1,
    strokecolor = :black,
    normalization = :pdf,
)
save(
    "figures/probability-escaping-trajectories.png",
    fig;
    resolution = (580, 300),
    px_per_unit = 3,
)
