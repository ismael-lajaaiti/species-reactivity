using CairoMakie
using DataFrames
using LinearAlgebra

A(n, sigma) = Float64.([i == j ? -1 : -abs(rand(Normal(0, sigma))) for i in 1:n, j in 1:n])
D(n, eps) = Diagonal(vcat([eps], fill(1, n - 1)))

sigma_values = [0.05, 0.1, 0.2]
n_range = 10:50
n_rep = 100
df = DataFrame(; n = Integer[], sigma = Float64[], determinant = Float64[])
for sigma in sigma_values, n in n_range, k in 1:n_rep
    push!(df, [n, sigma, abs(det(A(n, sigma)))])
end

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"n", ylabel = L"\det A")
colors = [:yellow, :blue, :purple]
strokewidth = 1
s_vec = []
for (sigma, color) in zip(sigma_values, colors)
    d = df[df.sigma.==sigma, :]
    d = combine(groupby(d, :n), :determinant => mean)
    s = scatter!(ax, d.n, d.determinant_mean; color, strokewidth)
    push!(s_vec, s)
end
Legend(
    fig[0, :],
    s_vec,
    [L"\sigma = %$sigma" for sigma in sigma_values];
    orientation = :horizontal,
)
save(
    "figures/determinant-random-matrix.png",
    fig;
    resolution = (450, 350),
    px_per_unit = 3,
)
