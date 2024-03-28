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
    collectivity(community::Community)

Compute the collectivity of the community.
The collectivity of the community is defined as the spectral radius of the interaction matrix.

For more information see: [Zelnick et al. (2024)](https://doi.org/10.1111/ele.14358).
"""
function collectivity(community::Community)
    maximum(norm.(eigvals(community.A)))
end

"""
    net_effects(community::Community)

Sum of the net effects taking place in the community.
Net effects are defined as the elements of the matrix `inv(A)`.
Returns the norm of the inverse of the interaction matrix.
"""
function net_effects(community::Community)
    norm(inv(community.A))
end

"""
    reactivity(com::Community)

Reactivity of the community.
"""
reactivity(com::Community) = reactivity(jacobian(com))
reactivity(j) = eigvals((j + transpose(j)) / 2)[end]

"""
    stability(com::Community)

Maximum real part of the eigenvalues of the Jacobian of the community.
"""
stability(com::Community) = stability(jacobian(com))
stability(j) = maximum(real(eigvals(j)))

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
    get_impact(c::Community, P::AbstractMatrix, f::Function)

Compute the impact of a perturbation `P` on a community `c` using the measure `f`.
Note that the measure `f` should be a function that takes a community as input and returns a scalar,
or a vector.
"""
function get_impact(c::Community, P::AbstractMatrix, f::Function)
    c_new = Community(c.A + P, c.r) # Community whose interactions have been perturbed.
    m = f(c) # Measure of interest for the original community.
    m_new = f(c_new) # Measure of interest for the perturbed community.
    norm(m - m_new) / norm(m)
end

function get_algebraic_impact(c::Community, P::AbstractMatrix, f::Function)
    c_new = Community(c.A + P, c.r) # Community whose interactions have been perturbed.
    m = f(c) # Measure of interest for the original community.
    m_new = f(c_new) # Measure of interest for the perturbed community.
    (m_new - m) / norm(m)
end

function synchrony(solution)
    u = solution.u
    u_matrix = reduce(hcat, u)
    total_species_var = (sum(std(u_matrix, dims=2)))^2
    total_community_var = var(sum(u_matrix, dims=1))
    total_community_var / total_species_var
end
