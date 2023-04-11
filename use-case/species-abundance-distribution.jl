using CairoMakie
using RareLotkaVolterra

S_values = [10, 30, 50]
std_values = [0.05, 0.1, 0.2]
n = 10 # Number of communities per richness and interaction standard deviation.
verbose = false

fig = Figure()

axes = Array{Any}(undef, length(S_values), length(std_values))
for (i, S) in enumerate(S_values), (j, std) in enumerate(std_values)
    @info "S = $S, std = $std"
    axes[i, j] = Axis(fig[i, j])
    j != 1 && hideydecorations!(axes[i, j]; grid = false)
    Label(
        fig[0, j],
        "σ = $std";
        tellwidth = false,
    )
    Label(
        fig[i, length(std_values)+1],
        "S = $S";
        rotation = -pi / 2,
        tellheight = false,
    )
    create_interaction_matrix(S) = random_competition_matrix(S, std)
    communities = create_communities(S, n; create_interaction_matrix, verbose)
    for c in communities[S]
        Neq = equilibrium_abundance(c)
        Neq_processed = log10.(sort(Neq / sum(Neq); rev = true))
        scatter!(Neq_processed; markersize = 8)
    end
end

linkyaxes!(reshape(axes, (length(S_values) * length(std_values),))...)

# Set x-axes label.
Label(
    fig[1+length(S_values), 1:length(std_values)],
    "Rank";
    tellwidth = false,
    width = 0.0,
)

# Set y-axes label.
Label(
    fig[2:end, 0],
    "Relative abundance (log)";
    rotation = pi / 2,
    tellheight = false,
    height = 0.0,
)

save(
    "figures/species-abundance-distribution.png",
    fig;
    resolution = (550, 550),
    px_per_unit = 3,
)
