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
    nonlinear_pb = ODEProblem(equilibrium_lotka_volterra, x0, tspan, [community, Neq])
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
