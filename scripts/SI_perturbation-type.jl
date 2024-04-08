using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using LinearAlgebra
using QuadGK
using Random
using SpeciesReactivity
using Statistics

Random.seed!(111) # For reproduciblity.

include("makie-theme.jl")

p_equal(S) = rand(Normal(), S)
p_proportional(S, yield) = p_equal(S) .* yield
p_inversely_proportional(S, yield) = p_equal(S) ./ yield

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n_communities = 1 # Number of communities.
mu = 0.0
sigma = 0.1
t_max = 1_000 # Duration of observation.
create_interaction_matrix(S) = random_competition_matrix(S, mu, sigma)
com = create_communities(S, n_communities; create_interaction_matrix)[S][1]
create_df() = DataFrame(;
    yield = Float64[],
    reactivity = Float64[],
    overall_response = Float64[],
    perturbation_type = String[],
)
n_perturbations = 100
intensity = 1
df = create_df()
A = com.A
yield = equilibrium_abundance(com)
reactivity = [get_reactivity(com.A, yield, i) for i in 1:S]
perturbation_dict = Dict(
    "Equal" => p_equal,
    "Proportional" => S -> p_proportional(S, yield),
    "Inversely proportional" => S -> p_inversely_proportional(S, yield),
)
for (p_type, p_func) in perturbation_dict
    @info "Perturbation type: $p_type"
    for k in 1:n_perturbations
        @info "Perturbation $k"
        x0 = p_func(S)
        z0 = x0 ./ yield
        x0 *= intensity / norm(z0)
        r = SpeciesReactivity.response(com, x0)
        for i in 1:S
            overall_response =
                quadgk(t -> abs(r.nonlinear(t; idxs = i)) / yield[i], 0, t_max)[1]
            # overall_response /= abs(x0[i] / yield[i]) # Normalization.
            push!(df, [yield[i], reactivity[i], overall_response, p_type])
        end
    end
end

processed_df = combine(
    groupby(df, [:reactivity, :yield, :perturbation_type]),
    :overall_response => mean => :overall_response,
    # :overall_response => maximum => :overall_response_worst,
    :overall_response =>
        (x -> std(x) / sqrt(n_perturbations)) => :overall_response_error,
)

with_theme(publication_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    c = fig[1, 3] = GridLayout()
    markersize = 8
    strokewidth = 0.5
    color = :grey
    p_type_ordered = ["Proportional", "Equal", "Inversely proportional"]
    p_type_formatted = ["Proportional", "Equal", "Inversely \n proportional"]
    for i in eachindex(p_type_ordered)
        p_type = p_type_ordered[i]
        title = p_type_formatted[i]
        sub_df = filter(row -> row.perturbation_type == p_type, processed_df)
        ax = Axis(
            fig[1, i];
            xlabel = "Species reactivity",
            ylabel = i == 1 ? "Species response intensity" : "",
            yscale = log10,
            title,
        )
        errorbars!(
            sub_df.reactivity,
            sub_df.overall_response,
            sub_df.overall_response_error;
            linewidth = 1,
            whiskerwidth = 3,
            color = :black,
        )
        scatter!(sub_df.reactivity, sub_df.overall_response; color, markersize, strokewidth)
    end
    for (layout, label) in zip([a, b, c], ["A", "B", "C"])
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    isdir("figures") || mkdir("figures")
    width = full_page_width * cm_to_pt
    height = width * 0.7 / width_height_ratio
    save_figure(
        "figures/SI_perturbation-type",
        # "/tmp/plot",
        fig,
        (width, height),
    )
end
