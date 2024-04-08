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

function random_matrix(S, mu, sigma)
    A = rand(Normal(mu, sigma), S, S)
    A - Diagonal(A) - I
end

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n_communities = 1 # Number of communities.
sigma = 0.1
t_max = 1_000 # Duration of observation.
n_perturbations = 100
create_df() = DataFrame(;
    yield = Float64[],
    reactivity = Float64[],
    overall_response = Float64[],
    mu = Float64[],
)
com_dict = Dict()
df = create_df()
mu_values = [0.02, -0.1, -0.05]
for mu in mu_values
    @info "Mean: $mu"
    create_interaction_matrix(S) = random_matrix(S, mu, sigma)
    com = create_communities(S, n_communities; create_interaction_matrix)[S][1]
    intensity = 0.6
    A = com.A
    com_dict[mu] = com
    yield = equilibrium_abundance(com)
    reactivity = [get_reactivity(com.A, yield, i) for i in 1:S]
    for k in 1:n_perturbations
        @info "Perturbation $k"
        x0 = proportional_perturbation(yield, intensity * sqrt(S), 1, true)
        r = SpeciesReactivity.response(com, x0)
        for i in 1:S
            overall_response =
                quadgk(t -> abs(r.nonlinear(t; idxs = i)) / yield[i], 0, t_max)[1]
            push!(df, [yield[i], reactivity[i], overall_response, mu])
        end
    end
end

processed_df = combine(
    groupby(df, [:reactivity, :yield, :mu]),
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
    d = fig[2, 1] = GridLayout()
    e = fig[2, 2] = GridLayout()
    f = fig[2, 3] = GridLayout()
    h = fig[3, 1] = GridLayout()
    i = fig[3, 2] = GridLayout()
    j = fig[3, 3] = GridLayout()
    markersize = 8
    strokewidth = 0.5
    color = :grey
    for (i, mu) in enumerate(mu_values)
        sub_df = filter(row -> row.mu == mu, processed_df)
        com = com_dict[mu]
        A = com.A
        frac_positive = count(>(0), A - Diagonal(A)) / (S * (S - 1))
        frac_positive = round(frac_positive; digits = 2)
        ax = Axis(
            fig[3, i];
            xlabel = "Species reactivity",
            ylabel = i == 1 ? "Species response \n intensity" : "",
            yscale = log10,
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
        ax = Axis(
            fig[2, i];
            xlabel = "Species relative yield",
            ylabel = i == 1 ? "Species reactivity" : "",
        )
        eta_min = minimum(sub_df.yield)
        eta_max = maximum(sub_df.yield)
        eta = range(eta_min, eta_max, 100)
        exp_reactivity_full = sqrt.(expected_reactivity_squared.(eta, Ref(com)))
        scatter!(sub_df.yield, sub_df.reactivity; color, markersize, strokewidth)
        lines!(
            eta,
            exp_reactivity_full;
            label = "Expectation",
            linewidth = 2,
            alpha = 0.7,
            color = palette[5],
        )
        ax = Axis(
            fig[1, i];
            xlabel = "Interaction strength",
            ylabel = i == 1 ? "Count" : "",
            title = "Fraction \n positive = $frac_positive",
        )
        hist!(
            vec(A);
            # bins = 20,
            color = :grey,
            strokewidth = 1,
            orientation = :vertical,
        )
        xlims!(-0.5, 0.3)
    end
    for (layout, label) in
        zip([a, b, c, d, e, f, h, i, j], ["A", "B", "C", "D", "E", "F", "H", "I", "J"])
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
    height = 2.5 * width * 0.7 / width_height_ratio
    save_figure(
        "figures/SI_positive-interaction",
        # "/tmp/plot",
        fig,
        (width, height),
    )
end
