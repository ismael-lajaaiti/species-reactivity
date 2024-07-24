using CairoMakie
using ColorSchemes
using DataFrames
using Distributions
using GLM
using LinearAlgebra
using QuadGK
using Random
using SpeciesReactivity
using Statistics

Random.seed!(14)

include("makie-theme.jl") # For reproducibility.

S = 50 # Community richness.
n = 50 # Number of communities.
sigma = 0.1
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
communities = create_communities(S, n; create_interaction_matrix)[S]
com = communities[2] # Get a community.
A = com.A
eta = equilibrium_abundance(com)
reactivity = [get_reactivity(A, eta, i) for i in 1:S]

function draw_carrying_capacities(eta)
    K1 = (1 .+ rand(Normal(0, 0.1), S)) ./ eta .^ 1.5
    K2 = (1 .+ rand(Normal(0, 0.1), S)) ./ eta
    K3 = (1 .+ rand(Normal(0, 0.1), S)) .* eta .^ 1
    [K1, K2, K3]
end
K_example = draw_carrying_capacities(eta)

n_perturbations = 100
intensity = 0.9
t_max = 1_000
df_list = Any[undef for _ in 1:n]
Threads.@threads for i in 1:n
    df_tmp = DataFrame(; K_index = Int64[], community_response = Float64[])
    c = communities[i]
    eta = equilibrium_abundance(c)
    K_vec = draw_carrying_capacities(eta)
    A = c.A
    for k in 1:n_perturbations
        @info "Perturbation $k"
        x0 = prop_perturbation(eta, intensity; no_extinction = true)
        r = SpeciesReactivity.response(c, x0)
        for (K_idx, K) in enumerate(K_example)
            delta_total_abundance(delta_eta) = abs(sum(K .* delta_eta))
            a = delta_total_abundance(x0) * t_max
            com_response =
                quadgk(t -> (delta_total_abundance(r.nonlinear(t)) / a), 0, t_max)[1]
            push!(df_tmp, [K_idx, com_response])
        end
    end
    df_list[i] = df_tmp
end
df = vcat(df_list...)

# Main figure.
with_theme(publication_theme) do
    fig = Figure()
    g = fig[1:3, 1:6] = GridLayout()
    g1 = g[1, 1:6] = GridLayout()
    g23 = g[2, 1:6] = GridLayout()
    g4 = fig[3, 1:6] = GridLayout()
    strokewidth = 0.5
    markersize = 9
    color = :grey
    # Panel A: reactivity vs. relative yield.
    ax = Axis(
        g[1, 1:6];
        aspect = AxisAspect(1.6),
        xlabel = "Species relative yield",
        ylabel = "Species reactivity",
    )
    scatter!(eta, reactivity; color, strokewidth, markersize)
    eta_linrange = LinRange(0, maximum(eta), 100)
    exp_reactivity_full = sqrt.(expected_reactivity_squared.(eta_linrange, Ref(com)))
    lines!(
        eta_linrange,
        exp_reactivity_full;
        color = palette[5],
        label = "Prediction",
        alpha = 0.7,
        linewidth = 2,
    )
    axislegend(; position = :rt)
    # Panel B and C: reactivity vs. abundance.
    for (i, color, K) in zip(1:3, palette[[2, 1, 3]], reverse(K_example))
        abundance = eta .* K
        relative_abundance = abundance / sum(abundance)
        # Reactivity vs. abundance: scatter plot.
        ylabel = i == 1 ? "Species reactivity" : ""
        i_ = i - 1
        col_idx = (2i_+1):(2i_+2)
        ax = Axis(
            g23[2, i];
            xlabel = "Species abundance",
            ylabel,
            titlesize = 12,
            tellwidth = true,
        )
        scatter!(
            relative_abundance,
            reactivity;
            color,
            strokewidth,
            markersize = 7,
            alpha = 0.6,
        )
        ax.xticks =
            round.(
                [
                    minimum(relative_abundance),
                    mean(extrema(relative_abundance)),
                    maximum(relative_abundance),
                ],
                digits = 3,
            )
    end
    Label(g23[1, 1], "Positive \n selection effect"; font = :bold)
    Label(g23[1, 2:3], "Negative selection effect"; font = :bold)
    # Panel C: histogram of overall community response. 
    ax4 = Axis(
        g[3, 1:6];
        aspect = AxisAspect(1.6),
        xlabel = "Log community overall \n response intensity",
        ylabel = "Density",
        tellwidth = false,
    )
    hist_vec = []
    for (i, color, K_idx) in zip(1:3, palette[[2, 1, 3]], [3, 2, 1])
        K = K_example[K_idx]
        tmp_hist = density!(
            ax4,
            log10.(df[df.K_index.==K_idx, :community_response]);
            color = (color, 0.3),
            strokecolor = color,
            strokewidth = 1.4,
        )
        # tmp_hist = hist!(
        #     ax4,
        #     # df[df.K_index.==K_idx, :community_response];
        #     log10.(df[df.K_index.==K_idx, :community_response]);
        #     normalization = :probability,
        #     color = (color, 0.6),
        #     bins = -3.2:0.1:0.1,
        #     strokewidth = 1,
        #     strokecolor = :white,
        #     label = i == 1 ? "Positive" : "Negative",
        # )
        push!(hist_vec, tmp_hist)
    end
    xlims!(-3.3, 0)
    # axislegend("Selection effect"; position = :rt, rowgap = 0)
    # Label panels.
    panels = [fig[1, 2], g23[2, 1], g23[2, 2], g23[2, 3], g[3, 2]]
    for (layout, label) in zip(panels, letters)
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    # xlims!(-0.1, 0.2)
    rowgap!(g23, 5)
    isdir("figures") || mkdir("figures")
    width = full_page_width * cm_to_pt
    height = width * (2 / 1.62)
    save_figure(
        "figures/02_selection-effect",
        # "/tmp/plot",
        fig,
        (width, height),
    )
end
