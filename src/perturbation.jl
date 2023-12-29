"""
    response(community::Community, x0; tspan = (0, 10_000))

Response of the `community` to an initial perturbation `x0` for a given time span `tspan`.

The response is quantified by the temporal evolution of the pertubation `x`.
Return the response of Lotka-Volterra system and
the reponse given by the Jacobian dynamics (i.e. the linearized dynamics).

```jldoctest
julia> community = create_communities(10, 1)[10][1] # Community of 10 species.
       x0 = 0.01 * rand(10) # Initial perturbation. x = response(community, x0) # Simulate recovery trajectory.
       x = response(community, x0; tspan = (0, 10_000))
       x.nonlinear(0) ≈ x.linear(0) ≈ x0
true

julia> all(abs.(x.nonlinear(1_000)) .< 1e-3) && all(abs.(x.linear(1_000)) .< 1e-3)
true
```

See also [`create_communities`](@ref).
"""
function response(community::Community, x0; tspan = (0, 10_000))
    Neq = equilibrium_abundance(community)
    N0 = Neq + x0
    @assert all(N0 .> 0) # Check that there is no extinction at t = 0.
    linear_pb = ODEProblem(jacobian_dynamics, x0, tspan, jacobian(community))
    nonlinear_pb = ODEProblem(equilibrium_lotka_volterra, x0, tspan, (community, Neq))
    (; linear = solve(linear_pb; alg = RK4()), nonlinear = solve(nonlinear_pb; alg = RK4()))
end

function isotrope_perturbation(Neq, intensity; no_extinction = false)
    S = length(Neq)
    function perturbation()
        x = rand(Normal(), S)
        (intensity / norm(x)) * x
    end
    no_extinction ? perturbation_no_extinction(perturbation, Neq) : perturbation()
end


function proportional_perturbation(
    Neq,
    expected_intensity,
    scaling_factor,
    no_extinction = false,
)
    function perturbation()
        x = isotrope_perturbation(Neq, 1)
        x = Diagonal(Neq) * x
        (expected_intensity / scaling_factor) * x
    end
    no_extinction ? perturbation_no_extinction(perturbation, Neq) : perturbation()
end

function expected_proportional_intensity(Neq; n = 100_000)
    mean([norm(proportional_perturbation(Neq, 1, 1)) for _ in 1:n])
end

function perturbation_no_extinction(create_perturbation, Neq; iter_max = 10^5)
    x = create_perturbation()
    iter = 0
    while any(Neq + x .<= 0) && iter < iter_max
        x = create_perturbation()
        iter += 1
    end
    iter == iter_max && throw(ErrorException("Maximum number of iteration reached."))
    x
end

amplitude(response, t_begin, t_end) = quadgk(t -> norm(response(t))^2, t_begin, t_end)[1]

function trajectory_error(x_true, x_hat; tspan = (0, 10_000))
    e = quadgk(t -> norm(x_true(t) - x_hat(t)), tspan[1], tspan[2]; rtol = 1e-2)[1]
    normalization = quadgk(t -> norm(x_true(t)), tspan[1], tspan[2]; rtol = 1e-2)[1]
    # normalization = norm(x_true(tspan[1]))
    # normalization = 1
    e / normalization
end

function transitional_regime_duration(
    x;
    threshold = 1e-2,
    t_init = 10_000,
    t_max = 1e5,
    delta_t = 1,
)
    t_end = t_init
    x0_norm = abs(x(0)) # Strength of initial perturbation.
    while t_end < t_max && abs(x(t_end)) > (threshold * x0_norm)
        t_end *= 10
    end
    @assert t_end <= t_max
    t_begin = t_end == t_init ? 0 : (t_end / 10)
    time_steps = values(t_begin:delta_t:t_end)
    x_sampled = abs.(x.(time_steps))
    x_sampled_reversed = reverse(x_sampled)
    for (i, x_t) in enumerate(x_sampled_reversed)
        x_t > (threshold * x0_norm) && return reverse(time_steps)[i]
    end
    return nothing
end

function nonlinear_excursion_data(x, N; tspan = (0, 100_000), delta_t = 1)
    time_steps = values(tspan[1]:delta_t:tspan[2])
    z(t) = x(t) / N
    z_sampled = abs.(z.(time_steps))
    z0 = abs(z(0))
    max_dist = maximum(z_sampled) - z0
    time = count(>=(z0), abs.(z_sampled)) * delta_t
    (; max_dist, time)
end

function integrate(fun, tspan; n_timesteps)
    trajectory_matrix = reduce(hcat, solution)
    S, n_timesteps = size(trajectory_matrix)
    integrals = zeros(S)
    for n_t in 1:(n_timesteps-1)
        delta_t = solution.t[n_t+1] - solution.t[n_t]
        for i in 1:S
            u_mean = (trajectory_matrix[n_t] + trajectory_matrix[n_t+1]) / 2
            integrals[i] += delta_t * u_mean
        end
    end
    integrals
end
