n_communities = 50 # Number of communities.
S = 30
sigma = 0.2
create_interaction_matrix(S) = random_competition_matrix(S, sigma)
communities = create_communities(S, n_communities; create_interaction_matrix)[S]
with_theme(p_theme) do
    fig = Figure()
    di_vec = Float64[]
    Di_vec = Float64[]
    di_prime_vec = Float64[]
    for com in communities
        eta = equilibrium_abundance(com)
        d_i = get_di(eta, 0)
        d_i_prime = get_di_prime(eta, 0)
        D_i = get_Di(eta, 0)
        push!(di_vec, d_i)
        push!(Di_vec, D_i)
        push!(di_prime_vec, d_i_prime)
    end
    println(di_vec)
    ax1 = Axis(fig[1, 1])
    scatter!(di_vec, Di_vec)
    ax2 = Axis(fig[1, 2])
    scatter!(di_prime_vec, Di_vec)
    save(
        # "figures/reactivity-yield-predictability.png",
        "/tmp/plot.png",
        fig;
        resolution = (600, 320),
        px_per_unit = 3,
    )
end
