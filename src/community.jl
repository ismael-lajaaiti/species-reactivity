"""
    Community(A::AbstractMatrix, r::AbstractVector)

Lotka-Volterra `Community` defined by its interaction matrix `A`
and its growth rate vector `r`.

# Example

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
    Community()

Create an empty [`Community`](@ref).

# Example

```jldoctest
julia> c = Community() # Create an empty community.
       (size(c.A), length(c.r)) == ((0, 0), 0)
true
```
"""
Community() = Community(Array{Any}(undef, 0, 0), [])


"""
    richness(community::Community)

Species richness of the [`Community`](@ref).

# Example

```jldoctest
julia> richness(Community()) == 0
true
```

See also [`Community`](@ref).
"""
richness(community::Community) = length(community.r)

"""
    equilibrium_abundance(community::Community)

Compute the abundance at equilibrium assuming a Lotka-Volterra system.

See also [`Community`](@ref).
"""
equilibrium_abundance(community::Community) = -inv(community.A) * fill(1, richness(community))

"""
    jacobian(community::Community)

Compute the Jacobian of the Lotka-Volterra system
given an interaction matrix `A` and the vector of abundances at equilibrium `Neq`.

See also [`Community`](@ref).
"""
function jacobian(community::Community)
    Neq = equilibrium_abundance(community)
    Diagonal(Neq .* community.r) * community.A
end

reactivity(j) = eigvals((j + transpose(j)) / 2)[end]

function nonlinear_reactivity(community::Community)
    Neq = equilibrium_abundance(community)
    j = community.A * Diagonal(Neq)
    eigvals((j + transpose(j)) / 2)[end]
end

"""
    is_stable(community::Community)

Is the `community` stable?

The community is said to be stable if the dominant eigenvalue
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
function lotka_volterra(N, community::Community, _)
    Diagonal(N .* community.r) * (1 .+ community.A * N)
end

"""
    equilibrium_lotka_volterra(x, p, _)

Lotka-Volterra differential equations of the species abundances `N` relative to the
equilibrium.
Thus, if the equilibrium is stable we expect that ||x|| tends to zero.
`p` is a collection whose first item contains the [`Community`](@ref)
and second item the species [`equilibrium_abundance`](@ref).
"""
function equilibrium_lotka_volterra(x, p, _)
    community, Neq = p
    Diagonal(community.r .* (Neq + x)) * (community.A * x)
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

function random_competition_matrix(S, mu, std)
    A = -abs.(rand(Normal(mu, std), S, S))
    A - Diagonal(A) - I
end

"""
    empty!(c::Community)

Remove all species of a [`Community`](@ref).

# Example

```jldoctest
julia> c = Community([0 -1; 0 -1], [1, 1])
       richness(c)
2

julia> empty!(c)
       richness(c)
0
```
"""
function Base.empty!(community::Community)
    empty_community = Community()
    community.A = empty_community.A
    community.r = empty_community.r
    return nothing
end

"""
    keep_surviving!(c::Community, surviving_species)

Remove species from the [`Community`](@ref)
that are not in the list of `surviving_species`.

# Example

```jldoctest
julia> c = Community([-1 0; 0 -2], [1, 2])
       surviving_species = []
       keep_surviving!(c, surviving_species)
       richness(c)
0

julia> c = Community([-1 0; 0 -2], [1, 2])
       surviving_species = [1]
       keep_surviving!(c, surviving_species)
       (c.A, c.r) == ([-1;;], [1])
true
```
"""
function keep_surviving!(c::Community, surviving_species)
    c.A = c.A[surviving_species, surviving_species]
    c.r = c.r[surviving_species]
    return nothing
end

"""
    assemble!(
        community::Community;
        tspan = (0, 1_000),
        extinction_threshold = 1e-6,
    )

Assemble the [`Community`](@ref), i.e. run the dynamics and filter extinct species.

## Keyword arguments

  - `tspan = (0, 1_000)`: duration of the simulation
  - `extinction_threshold = 1e-6`: species with lower abundance than this threshold
    at the end of a simulation are considered to be extinct

## Example

```jldoctest
julia> S = 20 # Species richness.
       sd = 0.1 # Interaction standard deviation.
       A = random_competition_matrix(S, sd)
       r = fill(1, S)
       c = Community(A, r)
       assemble!(c)
       richness(c) <= S
true
```

See also [`Community`](@ref) and [`keep_surviving!`](@ref).
"""
function assemble!(community::Community; tspan = (0, 1_000), extinction_threshold = 1e-6, return_surviving = false)
    S = richness(community)
    N0 = rand(Uniform(0.5, 1), S)
    problem = ODEProblem(lotka_volterra, N0, tspan, community)
    sol = solve(problem)
    if sol.retcode == :DtLessThanMin # Error during simulation, return empty community.
        empty!(community)
        return nothing
    end
    surviving_species = findall(>(extinction_threshold), sol.u[end])
    keep_surviving!(community, surviving_species)
    richness(community) == 0 && return nothing # No surviving species.
    Neq_surviving = equilibrium_abundance(community)
    if any(Neq_surviving .< 0) || !is_stable(community) # Check equilibrium.
        empty!(community)
    end
    return_surviving && return surviving_species
    return nothing
end

"""
    create_communities(
        Smin,
        Smax,
        n;
        create_interaction_matrix = S -> random_competition_matrix(S, 0.2),
        create_growth_rates = S -> fill(1, S),
        tspan = (0, 1_000),
        max_iter = 10^5,
    )

Create `n` Lotka-Volterra [`Community`](@ref) for each richness between `Smin` and `Smax`.

# Keyword arguments

  - `create_interaction_matrix = S -> [`random_competition_matrix`](@ref)(S, 0.2)`:
    function to create the interaction matrix given a richness
  - `create_growth_rates = S -> fill(1, S)`: function create the vector of
    growth rates
  - `tspan = (0, 1_000)` = time span of the simulations
  - `max_iter = 10^5` = maximum number of iteration before leaving the `while` loop
    if the communities have not been all created yet

# Example

```jldoctest
julia> c = create_communities(5, 10, 2) # Create 2 communities of S = 5, ..., 10 species.
       c8 = c[8][1] # First community of 8 species.
       richness(c8) == 8
true
```

See also [`Community`](@ref) and [`assemble!`](@ref).
"""
function create_communities(
    Smin,
    Smax,
    n;
    create_interaction_matrix = S -> random_competition_matrix(S, 0.2),
    create_growth_rates = S -> fill(1, S),
    tspan = (0, 1_000),
    max_iter = 10^5,
    verbose = false,
)
    # Communities start with richness in [Smin_pool, Smax_pool], with
    # Smin <= Smin_pool and Smax <= Smax_pool, as community richness can only
    # decrease during the assembly.
    Smax_pool = round(Integer, 1.5 * Smax)
    Smin_pool = Smin
    Srange = Smin:Smax
    community_dict = Dict([S => Community[] for S in Srange])
    is_filled(vec) = length(vec) == n
    all_filled(dict) = all(is_filled.(values(dict)))
    iter = 0
    while !all_filled(community_dict) && iter < max_iter
        for S in Smin_pool:Smax_pool
            r = create_growth_rates(S)
            A = create_interaction_matrix(S)
            community = Community(A, r)
            assemble!(community; tspan)
            Sfinal = richness(community)
            if Smin <= Sfinal <= Smax && !is_filled(community_dict[Sfinal])
                push!(community_dict[Sfinal], community)
            end
        end
        # Increase Smin_pool if all community of lower richness are already filled.
        idx = findfirst(S -> !is_filled(community_dict[S]), Smin:Smax)
        Smin_pool = isnothing(idx) ? nothing : Srange[idx]
        verbose && @info "Iteration $iter: $(length.(values(community_dict)))"
        iter += 1
    end
    community_dict
end

"""
    create_communities(S, n; kwargs...)

Create `n` Lotka-Volterra [`Community`](@ref) of size `S`.
"""
create_communities(S, n; kwargs...) = create_communities(S, S, n; kwargs...)

complexity(c::Community) = complexity(c.A)
function complexity(A::AbstractMatrix)
    norm(A - Diagonal(A))
end

departure_from_normality(c::Community) = departure_from_normality(c.A)
function departure_from_normality(A)
    sqrt(norm(A)^2 - sum(norm.(eigvals(A)).^2))
end

function shannon_diversity(Neq)
    p = Neq / sum(Neq)
    exp(- sum(p .* log.(p)))
end
