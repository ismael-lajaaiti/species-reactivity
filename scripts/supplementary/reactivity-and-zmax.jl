using DataFrames
using Distributions
using LinearAlgebra
using QuadGK
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
n_perturbation = 100
perturbation_intensity = 5
sigma_values = collect(0.05:0.05:0.25)
df = DataFrame(;
    sigma = Float64[], # Equivalent to community ID, as there is one community per sigma.
    species_id = Int64[],
    abundance = Float64[],
    max_nonlinear = Float64[],
)
A_dict = Dict()
Neq_dict = Dict()
for sigma in sigma_values
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    com = create_communities(S, n; create_interaction_matrix)[S][1] # Get a community.
    A = com.A
    A_dict[sigma] = A
    Neq = equilibrium_abundance(com)
    Neq_dict[sigma] = Neq
    for k in 1:n_perturbation
        @info "Perturbation $k"
        x0 = proportional_perturbation(Neq, perturbation_intensity, 1, true)
        r = response(com, x0)
        for i in 1:S
            nl(t) = abs(r.nonlinear(t; idxs = i)) / Neq[i]
            max_nl = find_max(nl; tspan = (0, 1_000))
            push!(df, [sigma, i, Neq[i], max_nl])
        end
    end
end
# df = df[df.return_time.<10_000, :]
pdf = combine(
    groupby(df, [:sigma, :species_id, :abundance]),
    :max_nonlinear => mean => :max_nonlinear_mean,
)
remove_diagonal(A) = A - Diagonal(A)
transform!(
    pdf,
    [:sigma, :species_id] =>
        ByRow(
            (sigma, i) ->
                sqrt(sum((remove_diagonal(A_dict[sigma])[i, :] .* Neq_dict[sigma]) .^ 2)),
        ) => :reactivity,
    [:sigma, :species_id] =>
        ByRow(
            (sigma, i) ->
                sqrt(
                    Neq_dict[sigma][i]^2 * (pi / 2 - 1) +
                    sum((remove_diagonal(A_dict[sigma])[i, :] .* Neq_dict[sigma]) .^ 2),
                ) - Neq_dict[sigma][i],
        ) => :reactivity_bis,
)
pdf_per_communities = combine(
    groupby(pdf, :sigma),
    :max_nonlinear_mean => mean => :max_nonlinear_mean,
    :reactivity => mean,
)
pdf_per_species = pdf[pdf.sigma.==sigma_values[end], :]

with_theme(p_theme) do
    fig = Figure()
    a = fig[1, 1] = GridLayout()
    b = fig[1, 2] = GridLayout()
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"Reactivity, $<R_{0, i}>_i$",
        ylabel = L"Maximum nonlinearity, $<\mathcal{z}_{\mathrm{max}, i}>_i$",
        title = "Community level",
    )
    scatter!(
        pdf_per_communities.reactivity_mean,
        pdf_per_communities.max_nonlinear_mean;
        color = :black,
    )
    ax2 = Axis(
        fig[1, 2];
        xlabel = L"Reactivity, $R_{0, i}$",
        ylabel = L"Maximum nonlinearity, $\mathcal{z}_{\mathrm{max}, i}$",
        title = "Species level",
    )
    scatter!(
        pdf_per_species.reactivity,
        pdf_per_species.max_nonlinear_mean;
        color = pdf_per_species.abundance,
    )
    colorrange = extrema(pdf_per_species.abundance)
    Colorbar(fig[1, 3]; label = L"Species abundance, $\eta_i$", colorrange)
    save(
        "figures/reactivity-and-abundance.png",
        fig;
        resolution = (620, 320),
        px_per_unit = 3,
    )
end
