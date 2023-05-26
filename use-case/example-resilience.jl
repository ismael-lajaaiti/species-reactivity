using CairoMakie
using DataFrames
using LaTeXStrings
using LinearAlgebra
using RareLotkaVolterra
using Statistics
import ColorSchemes: leonardo

S = 50 # Community richness.
n = 1
sigma = 0.1 # Interspecific interaction strengths.
verbose = false
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
c = create_communities(S, n; create_interaction_matrix, verbose)[S][1]
Neq = equilibrium_abundance(c)
x0 = minimum(Neq) / 2 * ones(S)
r = response(c, x0)
timesteps = LinRange(0, 200, 10_000) #= 10 .^ LinRange(-3, 4, 100) =#

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Time", ylabel = "Distance to equilibrium (log)")
trajectories = reduce(hcat, r.nonlinear.(timesteps))
l_rare = lines!(
    timesteps,
    log10.(norm.(trajectories[argmin(Neq), :]) / abs(x0[begin]));
    color = :blue,
)
l_com = lines!(
    timesteps,
    log10.(vec(mean(abs.(trajectories); dims = 1)) / abs(x0[begin]));
    color = :red,
)
text!(L"- R_\infty t"; position = (150, -0.1))
text!(L"- R_\infty t"; position = (150, -1.2))
Legend(
    fig[0, :],
    [l_com, l_rare],
    ["Community (mean)", "Rarest species"];
    orientation = :horizontal,
)
save("figures/example-resilience2.png", fig; resolution = (450, 350), px_per_unit = 3)

rare = argmin(Neq) # Index of the least abundant species.
j = jacobian(c)


fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"N_\mathrm{min}", ylabel = L"\lambda_\mathrm{dom}")
# ax2 = Axis(fig[1, 2]; xlabel = L"N_\mathrm{min}", ylabel = L"\mathbf{v}_\mathrm{dom}_{\mathbf{N}^*}")
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
