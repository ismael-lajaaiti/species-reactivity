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

Random.seed!(112)

include("makie-theme.jl") # For reproducibility.

S = 30
n = 1
sigma = 0.25
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
eta = equilibrium_abundance(com)

# Compute the probability of generating an isotrope perturbation outside the linearity
# domain for different `intensity_values`.
S = 50 # Community richness.
n = 1 # Number of communities.
sigma = 0.15
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
A = com.A
eta = equilibrium_abundance(com)
reactivity = [get_reactivity(A, eta, i) for i in 1:S]

K1 = (1 .+ rand(Normal(0, 0.1), S)) ./ eta .^ 1.25
K2 = (1 .+ rand(Normal(0, 0.1), S)) ./ eta
K3 = (1 .+ rand(Normal(0, 0.1), S)) .* eta
K_vec = [K1, K2, K3]
n_perturbations = 1_000
intensity = 0.6
t_max = 1_000
df = DataFrame(; K_index = Int64[], community_response = Float64[])
for k in 1:n_perturbations
    @info "Perturbation $k"
    A = com.A
    x0 = proportional_perturbation(eta, intensity * sqrt(S), 1, true)
    r = SpeciesReactivity.response(com, x0)
    for (K_idx, K) in enumerate(K_vec)
        delta_total_abundance(delta_eta) = abs(sum(K .* delta_eta))
        a = delta_total_abundance(x0) * t_max
        com_response = quadgk(t -> (delta_total_abundance(r.nonlinear(t)) / a), 0, t_max)[1]
        push!(df, [K_idx, com_response])
    end
end

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
        label = "expectation",
        alpha = 0.7,
        linewidth = 2,
    )
    axislegend(; position = :rt)
    # Panel B and C: reactivity vs. abundance.
    for (i, color, K) in zip(1:3, palette[[2, 1, 3]], [K3, K2, K1])
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
            alpha = 0.8,
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
        ylabel = "Frequency",
        tellwidth = false,
    )
    hist_vec = []
    for (i, color, K_idx) in zip(1:3, palette[[2, 3, 1]], [3, 1, 2])
        K = K_vec[K_idx]
        tmp_hist = hist!(
            ax4,
            log10.(df[df.K_index.==K_idx, :community_response]);
            normalization = :probability,
            color = (color, 0.7),
            bins = -3.2:0.2:1,
            strokewidth = 1,
            strokecolor = :white,
            label = i == 1 ? "Positive" : "Negative",
        )
        push!(hist_vec, tmp_hist)
    end
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
    rowgap!(g23, 5)
    isdir("figures") || mkdir("figures")
    width = full_page_width * cm_to_pt
    height = width * (2 / 1.62)
    save(
       # "figures/02_selection-effect.pdf",
        "/tmp/plot.png",
        fig;
        size = (width, height),
        pt_per_unit = 1,
    )
end

# Supplementary figure.
with_theme(publication_theme) do
    fig = Figure()
    strokewidth = 0.5
    markersize = 9
    n_panels = 10
    panels = Any[nothing for _ in 1:n_panels]
    panels[1] = fig[1, 2] = GridLayout()
    panels[2] = fig[3, 1] = GridLayout()
    panels[3] = fig[3, 2] = GridLayout()
    panels[4] = fig[3, 3] = GridLayout()
    panels[5] = fig[4, 1] = GridLayout()
    panels[6] = fig[4, 2] = GridLayout()
    panels[7] = fig[4, 3] = GridLayout()
    panels[8] = fig[5, 1] = GridLayout()
    panels[9] = fig[5, 2] = GridLayout()
    panels[10] = fig[5, 3] = GridLayout()
    title_vec = [L"\mathrm{Cov}(N, η)<0", L"\mathrm{Cov}(N, η)=0", L"\mathrm{Cov}(N, η)>0"]
    color = :grey
    for (i, K) in enumerate(K_vec)
        title = title_vec[i]
        abundance = eta .* K
        relative_abundance = abundance / sum(abundance)
        # Abundance vs. relative yield.
        ylabel = i == 1 ? "Abundance (N)" : ""
        ax = Axis(fig[3, i]; xlabel = "Relative yield (η)", ylabel, title)
        scatter!(ax, eta, relative_abundance; color, strokewidth, markersize)
        # Reactivity vs. abundance: scatter plot.
        ylabel = i == 1 ? "Reactivity" : ""
        ax2 = Axis(fig[4, i]; xlabel = "Abundance", ylabel)
        scatter!(ax2, relative_abundance, reactivity; color, strokewidth, markersize)
        # Reactivity vs. abundance: linear fit.
        df_tmp = DataFrame(; reactivity, relative_abundance, eta)
        ols = lm(@formula(reactivity ~ relative_abundance), df_tmp)
        p_value = ftest(ols.model).pval
        a, b = coef(ols)
        xlims = [minimum(relative_abundance), maximum(relative_abundance)]
        ylims = a .+ b .* xlims
        linestyle = p_value > 0.05 ? :dot : :solid
        lines!(xlims, ylims; color = :skyblue, linestyle, linewidth = 2, label = "fit")
        i == 3 && axislegend(ax2; labelsize = 10)
        # Histogram of overall community responses.
        limits = extrema(log10.(df.community_response))
        ylabel = i == 1 ? "Frequency" : ""
        ax = Axis(fig[5, i]; xlabel = "Log community overall \n response intensity", ylabel)
        hist!(
            log10.(df[df.K_index.==i, :community_response]);
            normalization = :probability,
            color = :grey,
            # strokewidth = 1,
            bins = 15,
            strokewidth = 1,
            strokecolor = :white,
        )
        vlines!(
            mean(log10.(df[df.K_index.==i, :community_response]));
            color = :black,
            linewidth = 2,
            linestyle = :dash,
            label = "mean",
        )
        xlims!(limits...)
        i == 3 && axislegend(ax; labelsize = 10)
    end
    ax = Axis(fig[1, 2]; xlabel = "Relative yield", ylabel = "Reactivity")
    scatter!(eta, reactivity; color, strokewidth, markersize)
    eta_linrange = LinRange(0, maximum(eta), 100)
    exp_reactivity_full = sqrt.(expected_reactivity_squared.(eta_linrange, Ref(com)))
    lines!(
        eta_linrange,
        exp_reactivity_full;
        color = :orange,
        label = "expectation",
        alpha = 0.7,
        linewidth = 2,
    )
    axislegend(; position = :rt)
    Label(fig[2, 1:2], "Selection effect < 0"; tellwidth = false, fontsize = 15)
    Label(fig[2, 3], "Selection effect ≥ 0"; tellwidth = false, fontsize = 15)
    # Label panels.
    for (layout, label) in zip(panels, letters)
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
    isdir("figures") || mkdir("figures")
    width = 1.1 * full_page_width * cm_to_pt
    height = width * 2 / width_height_ratio
    save(
        "figures/SI_02_selection-effect.pdf",
        # "/tmp/plot.png",
        fig;
        size = (width, height),
        pt_per_unit = 1,
    )
end
