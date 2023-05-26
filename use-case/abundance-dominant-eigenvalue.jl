using CairoMakie
using DataFrames
using LaTeXStrings
using LinearAlgebra
using RareLotkaVolterra

S = 50 # Community richness.
n = 100 # Number of communities.
sigma_values = [0.1] # Interspecific interaction strengths.
verbose = false # Show @info while generating communities.

function mass(x, weight)
    x /= norm(x)
    weight /= sum(weight)
    (dot(x .^ 2, weight) - minimum(weight)) / (maximum(weight) - minimum(weight))
end

df = DataFrame(;
    sigma = Float64[],
    Neq_rare = Float64[], # Relative abundance of the rarest species.
    lambda_dom = Float64[], # Dominant eigenvalue of Jacocbian.
    pred_lambda_dom = Float64[], # Predicted dominant eigenvalue of Jacocbian.
    coef_rare_squared = Float64[], # Squared coefficient of the dominant eigenvector
    # of the rarest species.
    mass = Float64[],
)
for sigma in sigma_values
    create_interaction_matrix(S) = random_competition_matrix(S, sigma)
    communities = create_communities(S, n; create_interaction_matrix, verbose)
    for c in communities[S]
        Neq = equilibrium_abundance(c)
        rare = argmin(Neq) # Index of the least abundant species.
        j = jacobian(c)
        lambda_dom = real(eigvals(j)[end])
        all_but_rare = [i for i in 1:S if i != rare]
        A_without_rare = c.A[all_but_rare, all_but_rare]
        pred_lambda_dom = Neq[rare] * (det(c.A) / det(A_without_rare))
        vec_dom = eigvecs(j)[:, end]
        coef_rare_squared = norm(vec_dom[rare])^2
        m = mass(norm.(vec_dom), Neq)
        push!(df, [sigma, Neq[rare], lambda_dom, pred_lambda_dom, coef_rare_squared, m])
    end
    @info "sigma = $sigma done."
end

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"N_\mathrm{min}", ylabel = L"\lambda_\mathrm{dom}")
# ax2 = Axis(fig[1, 2]; xlabel = L"N_\mathrm{min}", ylabel = L"<\mathbf{v}_\mathrm{dom}>_{\mathbf{N}^*}")
ax3 = Axis(
    fig[1, 2];
    xlabel = L"\lambda_\mathrm{dom}",
    ylabel = L"\hat{\lambda}_\mathrm{dom}",
)
colors = [:yellow, :blue, :purple]
strokewidth = 1
s_vec = []
for (sigma, color) in zip(sigma_values, colors)
    d = df[df.sigma.==sigma, :]
    s = scatter!(ax, d.Neq_rare, d.lambda_dom; color, strokewidth)
    # scatter!(ax2, d.Neq_rare, d.mass; color, strokewidth)
    scatter!(ax3, d.lambda_dom, d.pred_lambda_dom; color, strokewidth)
    push!(s_vec, s)
end
# Legend(
#     fig[0, :],
#     s_vec,
#     [L"\sigma = %$sigma" for sigma in sigma_values];
#     orientation = :horizontal,
# )
save(
    "figures/abundance-dominant-eigenvalue-new.png",
    fig;
    resolution = (500, 250),
    px_per_unit = 3,
)

fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel = L"N_\mathrm{min}",
    ylabel = L"S|\mathbf{v}_\mathrm{dom}^\mathrm{rare}|^2",
)
colors = [:yellow, :blue, :purple]
strokewidth = 1
s_vec = []
for (sigma, color) in zip(sigma_values, colors)
    d = df[df.sigma.==sigma, :]
    s = scatter!(ax, d.Neq_rare, S .* d.coef_rare_squared; color, strokewidth)
    push!(s_vec, s)
end
# Legend(
#     fig[0, :],
#     s_vec,
#     [L"\sigma = %$sigma" for sigma in sigma_values];
#     orientation = :horizontal,
# )
save(
    "figures/abundance-dominant-eigenvector.png",
    fig;
    resolution = (300, 250),
    px_per_unit = 3,
)
