using CairoMakie
using ColorSchemes
using Distributions
using LaTeXStrings
using RareLotkaVolterra
using Statistics

S = 10
c = create_communities(S, 1)[S][1]
Neq = equilibrium_abundance(c)
x0_unscaled = abs.(rand(Normal(), S))
colors = colorschemes[:tol_light]
fig = Figure()
for (i, intensity) in enumerate([0.1, 1.0])
    ax = Axis(
        fig[1, i];
        xlabel = L"\log t",
        ylabel = L"\mathbf{x}",
        title = "intensity = $intensity",
    )
    x0 = intensity * mean(Neq) * x0_unscaled
    println(x0)
    x = response(c, x0)
    for (sp, col) in zip(1:S, colors)
        lines!(
            log10.(x.nonlinear.t[begin+1:end]),
            x.nonlinear[sp, begin+1:end] .- Neq[sp];
            color = col,
        )
        lines!(
            log10.(x.linear.t[begin+1:end]),
            x.linear[sp, begin+1:end];
            color = col,
            linestyle = :dash,
        )
    end
end
save(
    "/tmp/plot.png", # Do not forget to replace with the relevant path.
    fig;
    resolution = (550, 350),
    px_per_unit = 3,
)

# colors = [:yellow, :blue, :purple]
# strokewidth = 1
# s_vec = []
# for (sigma, color) in zip(sigma_values, colors)
#     d = df[df.sigma.==sigma, :]
#     s = scatter!(ax, d.p_rare, d.lambda_dom; color, strokewidth)
#     scatter!(ax2, d.p_rare, d.mass; color, strokewidth)
#     push!(s_vec, s)
# end
# Legend(
#     fig[0, :],
#     s_vec,
#     [L"\sigma = %$sigma" for sigma in sigma_values];
#     orientation = :horizontal,
# )
