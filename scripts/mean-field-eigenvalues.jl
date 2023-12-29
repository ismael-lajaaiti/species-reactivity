using CairoMakie
using DataFrames
using Distributions
using LinearAlgebra

A(n, sigma) = [i == j ? -1 : -abs(rand(Normal(0, sigma))) for i in 1:n, j in 1:n]
D(n, eps) = Diagonal(vcat([eps], fill(1, n - 1)))
function mass(x, weight)
    x /= norm(x)
    weight /= sum(weight)
    (dot(x .^ 2, weight) - minimum(weight)) / (maximum(weight) - minimum(weight))
end

df = DataFrame(;
    sigma = Float64[],
    epsilon = Float64[],
    dominant_eigval = Float64[],
    dominant_mass = Float64[],
)
n = 50
sigma_values = [0.05, 0.1, 0.2]
epsilon_values = LinRange(0.01, 0.5, 10)
n_rep = 100
for sigma in sigma_values, epsilon in epsilon_values, i in 1:n_rep
    M = D(n, epsilon) * A(n, sigma)
    M = Float64.(M) # Necessary convertion to compute eigen values.
    dominant_eigval = real(eigvals(M)[end])
    Neq = vcat([epsilon], fill(1, n - 1))
    dominant_mass = mass(norm.(eigvecs(M)[:, end]), Neq)
    push!(df, [sigma, epsilon, dominant_eigval, dominant_mass])
end

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"\epsilon", ylabel = L"\lambda_\mathrm{dom}")
ax2 = Axis(
    fig[1, 2];
    xlabel = L"\epsilon",
    ylabel = L"<\mathbf{v}_\mathrm{dom}>_{\mathbf{N}^*}",
)
colors = [:yellow, :blue, :purple]
strokewidth = 1
s_vec = []
for (sigma, color) in zip(sigma_values, colors)
    d = df[df.sigma.==sigma, :]
    d = combine(groupby(d, :epsilon), :dominant_eigval => mean, :dominant_mass => mean)
    s = scatter!(ax, d.epsilon, d.dominant_eigval_mean; color, strokewidth)
    scatter!(ax2, d.epsilon, d.dominant_mass_mean; color, strokewidth)
    push!(s_vec, s)
end
Legend(
    fig[0, :],
    s_vec,
    [L"\sigma = %$sigma" for sigma in sigma_values];
    orientation = :horizontal,
)
save(
    "figures/abundance-dominant-eigenvalue-simplified.png",
    fig;
    resolution = (500, 250),
    px_per_unit = 3,
)
