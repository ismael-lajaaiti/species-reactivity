"""
    Community(A::AbstractMatrix, r::AbstractVector)

Lotka-Volterra `Community` defined by its interaction matrix `A`
and its growth rate vector `r`.


```jldoctest
julia> A = [-1 0; 0 -1]
       r = [1, 1]
       Community(A, r)
Community([-1 0; 0 -1], [1, 1])
```

See also [`lotka_volterra_dynamics`](@ref), [`equilibrium_abundance`](@ref)
and [`jacobian`](@ref).
"""
mutable struct Community
    A::AbstractMatrix
    r::AbstractVector
    function Community(A, r)
        @assert size(A, 1) == size(A, 2) == length(r)
        new(A, r)
    end
end

"""
    empty_community()

Return an empty [`Community`](@ref).
"""
empty_community() = Community(Array{Any}(undef, 0, 0), [])

"""
    richness(community::Community)

Species richness of the [`Community`](@ref).
"""
richness(community::Community) = length(community.r)

"""
    isempty(community::Community)

Is the community empty, i.e. has a zero species richness?

See also [`richness`](@ref) and [`Community`](@ref).
"""
Base.isempty(community::Community) = richness(community) == 0

"""
    equilibrium_abundance(community::Community)

Compute the abundance at equilibrium assuming a Lotka-Volterra system.

See also [`Community`](@ref).
"""
equilibrium_abundance(community::Community) = -inv(community.A) * community.r

"""
    jacobian(community::Community)

Compute the Jacobian of the Lotka-Volterra system
given an interaction matrix `A` and the vector of abundances at equilibrium `Neq`.

See also [`Community`](@ref).
"""
function jacobian(community::Community)
    Neq = equilibrium_abundance(community)
    Diagonal(Neq) * community.A
end

"""
    is_stable(community::Community)

Is the `community` stable?

The community is said to be stable if the dominant eigen value
of its [`jacobian`](@ref) is negative.

See also [`Community`](@ref).
"""
function is_stable(community::Community)
    J = jacobian(community)
    maximum(real.(eigvals(J))) < 0
end

"""
    lotka_volterra_dynamics(N, community::Community, _)

Lotka-Volterra differential equations of the species abundances `N`.
"""
function lotka_volterra_dynamics(N, community::Community, _)
    Diagonal(N) * (community.r + community.A * N)
end

"""
    jacobian_dynamics(x, J, _)

Linearized dynamics given by the `J`acobian of a perturbation `x`.
"""
jacobian_dynamics(x, J, _) = J * x

"""
    random_competition_matrix(S, std)

Create randomly an interaction with competitive interactions
(i.e. negative coefficients).

The output matrix has a size of (`S`, `S`) where `S` is the species richness,
and the coefficients are taken in a normal law centered in 0
and of standard deviation `std`.
To ensure that coeffients are negative we apply the following transformation:
`x -> - abs(x)`.
The diagonal elements are set to `-1`.
"""
function random_competition_matrix(S, std)
    A = -abs.(rand(Normal(0, std), S, S))
    A - Diagonal(A) - I
end

"""
    assemble!(
        community::Community;
        tspan = (0, 100),
        extinction_threshold = 1e-4,
    )

Assemble the `community`.

Run the the dynamics for the defined `tspan` and once the equilibrium is reached
remove extinct species, i.e. species whose abundance is below the `extinction_threshold`.

See also [`Community`](@ref).
"""
function assemble!(
    community::Community;
    tspan = (0, 1_000),
    extinction_threshold = 1e-4,
)
    S = richness(community)
    N0 = rand(Uniform(0.1, 1), S) # initial abundances
    problem = ODEProblem(lotka_volterra_dynamics, N0, tspan, community)
    sol = solve(problem; reltol = 1e-8)
    # If error during the simulation, return prematurly an empty community.
    sol.retcode == :DtLessThanMin && return empty_community()
    is_extinct = sol[end] .< extinction_threshold # extinct species below threshold
    surviving_species = (1:S)[.!is_extinct]
    # If not surviving species, return prematurly nothing.
    isempty(surviving_species) && return empty_community()
    # Update community by keeping only surviving species and check that is stable.
    community.A = community.A[surviving_species, surviving_species]
    community.r = community.r[surviving_species]
    Neq_surviving = equilibrium_abundance(community)
    # If a negative biomass or unstable equilibrium return an empty community.
    if !all(Neq_surviving .> 0) || !is_stable(community)
        empty = empty_community()
        community.A = empty.A
        community.r = empty.r
    end
end

"""
    create_communities(
        Smin,
        Smax,
        n_communities;
        create_interaction_matrix = S -> random_competition_matrix(S, 0.1),
        create_growth_rates = S -> fill(1, S),
        tspan = (0, 1_000),
        max_iter = 10^5,
    )

Create `n_communities` for each richness between `Smin` and `Smax`.

## Keyword arguments

  - `create_interaction_matrix = S -> [`random_competition_matrix`](@ref)(S, 0.2)`:
    function that creates the interaction matrix given a richness.
  - `create_growth_rates = S -> fill(1, S)`: function that creates the vector of
    growth rates.
  - `tspan = (0, 1_000)` = time span of the simulations.
  - `max_iter = 10^5` = maximum number of iteration before quiting the while loop
    if the communities have not been all created yet.

See also [`Community`](@ref).
"""
function create_communities(
    Smin,
    Smax,
    n_communities;
    create_interaction_matrix = S -> random_competition_matrix(S, 0.2),
    create_growth_rates = S -> fill(1, S),
    tspan = (0, 1_000),
    max_iter = 10^5,
)
    # Initially the communities will start with a richness between Smin_pool and
    # Smax_pool, with Smax_pool >= Smax has we know that during the assembly
    # the community will loose species.
    # Thus by starting with an initial higher we increase our chances to have a final
    # diversity inside the target range [Smin; Smax].
    Smax_pool = round(Integer, 1.5 * Smax)
    Smin_pool = Smin
    Srange = Smin:Smax
    community_dict = Dict([S => Community[] for S in Srange])
    iter = 0
    is_filled(vec, n) = length(vec) == n
    all_community_filled(dict, n) = all(is_filled.(values(dict), n))
    while !all_community_filled(community_dict, n_communities) && iter < max_iter
        for S in Smin_pool:Smax_pool
            r = create_growth_rates(S)
            A = create_interaction_matrix(S)
            community = Community(A, r)
            assemble!(community; tspan = tspan)
            Sfinal = richness(community)
            if Smin <= Sfinal <= Smax && !is_filled(community_dict[Sfinal], n_communities)
                push!(community_dict[Sfinal], community)
            end
        end
        # Increase Smin_pool if all community of lower richness are already filled.
        Smin_pool = findfirst(S -> !is_filled(community_dict[S], n_communities), Smin:Smax)
        iter += 1
    end
    community_dict
end

function create_communities(S, n_communities; kwargs...)
    create_communities(S, S, n_communities; kwargs...)
end
