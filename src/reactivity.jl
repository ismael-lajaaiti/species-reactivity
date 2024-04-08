remove_diagonal(A) = A - Diagonal(A)
get_reactivity(A, Neq, i) = sqrt(sum((remove_diagonal(A)[i, :] .* Neq) .^ 2))
function get_reactivity(c::Community)
    S = richness(c)
    eta = equilibrium_abundance(c)
    A = remove_diagonal(c.A)
    reactivity(A, eta, i) = sqrt(sum((A[i, :] .* eta) .^ 2))
    [reactivity(A, eta, i) for i in 1:S]
end

get_Di(eta, eta_i) = (sum(eta) - eta_i)^2 / (norm(eta)^2 - eta_i^2)
get_di(eta, eta_i) = (norm(eta)^2 - eta_i^2)^2 / (norm(eta, 4)^4 - eta_i^4)
get_di_prime(eta, eta_i) = (norm(eta)^2 - eta_i^2)^3 / (norm(eta, 3)^3 - eta_i^3)^2

function expected_reactivity_squared(eta_i, com)
    S = richness(com)
    d = S - 1
    A = remove_diagonal(com.A)
    eta = equilibrium_abundance(com)
    mean_ai = sum(A) / (S * d)
    mean_ai2 = sum(A .^ 2) / (S * d)
    D_i = get_Di(eta, eta_i)
    d_i = get_di(eta, eta_i)
    d_i_prime = get_di_prime(eta, eta_i)
    a =
        mean_ai2 * (1 - 1 / d_i) -
        mean_ai^2 * (2 * sqrt(D_i / d_i_prime) - D_i / d_i - 1 / d_i) +
        2 *
        mean_ai *
        ((1 - eta_i) / (sum(eta) - eta_i)) *
        (D_i / d_i - sqrt(D_i / d_i_prime)) +
        (D_i / d_i) * ((1 - eta_i) / (sum(eta) - eta_i))^2
    a * (norm(eta)^2 - eta_i^2)
end

function expected_reactivity_squared(eta_i, com, C)
    S = richness(com)
    d = S - 1
    A = remove_diagonal(com.A)
    eta = equilibrium_abundance(com)
    mean_ai = C * mean(A[A.!=0])
    mean_ai2 = C * mean(A[A.!=0] .^ 2)
    D_i = get_Di(eta, eta_i)
    d_i = get_di(eta, eta_i)
    d_i_prime = get_di_prime(eta, eta_i)
    a =
        mean_ai2 * (1 - 1 / d_i) -
        mean_ai^2 * (2 * sqrt(D_i / d_i_prime) - D_i / d_i - 1 / d_i) +
        2 *
        mean_ai *
        ((1 - eta_i) / (sum(eta) - eta_i)) *
        (D_i / d_i - sqrt(D_i / d_i_prime)) +
        (D_i / d_i) * ((1 - eta_i) / (sum(eta) - eta_i))^2
    a * (norm(eta)^2 - eta_i^2)
end

function expected_reactivity_squared_naive(eta_i, com)
    S = richness(com)
    d = S - 1
    A = remove_diagonal(com.A)
    eta = equilibrium_abundance(com)
    mean_ai2 = sum(A .^ 2) / (S * d)
    mean_ai2 * (norm(eta)^2 - eta_i^2)
end

