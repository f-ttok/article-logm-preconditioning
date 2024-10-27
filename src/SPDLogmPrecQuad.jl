module SPDLogmPrecQuad
using LinearAlgebra
using SuiteSparse
using FastGaussQuadrature

export logm_gl, logm_pgl, logm_de, logm_pde

struct LogmGLResult
    X
    m
end

"""
Computing log(A)B using the Gauss--Legendre quadrature.
"""
function logm_gl(A, λ_max, λ_min; B=I, ϵ=1e-12, m::Union{Nothing,Int}=nothing)
    c = 1 / sqrt(λ_max * λ_min)
    κ = λ_max / λ_min
    if isnothing(m)
        m = select_nof_abscissas_gl(κ^(1 / 2), ϵ)
    end
    F(t) = ((1 + t) * c * A + (1 - t) * I) \ B
    t, w = gausslegendre(m)
    G = zero(B)
    for k in 1:m
        G += w[k] * F(t[k])
    end

    G = (c * A - I) * G - log(c) * B
    return LogmGLResult(G, m)
end

"""
Select the number of abscissas for the GL quadrature for log(λ)
(λ> 0) so that the error is less than ϵ.
"""
function select_nof_abscissas_gl(λ, ϵ)
    f(t) = (λ - 1) / ((1 + t) * λ + (1 - t))
    err = Inf
    exact = log(λ)
    m_max = 1000
    m_min = 2
    while m_max - m_min > 1
        m = (m_max + m_min) ÷ 2
        t, w = gausslegendre(m)
        approx = sum(@. w * f(t))
        err = abs(exact - approx)
        if err > ϵ
            m_min = m
        else
            m_max = m
        end
    end
    return m_max
end

"""
Computing log(A)B using the preconditioned GL quadrature.
"""
function logm_pgl(A, λ_max, λ_min; B=I, ϵ=1e-12, m::Union{Nothing,Int}=nothing)
    κ = λ_max / λ_min
    c = 1 / sqrt(λ_max * λ_min)
    cc = sqrt((sqrt(κ) + 1) * (1 / sqrt(κ) + 1)) # = c'
    F1(t) = (((cc - 1) * t + 1 + cc) * c * A + (1 - t) * I) \ B
    F2(t) = ((1 - t) * c * A + ((cc - 1) * t + cc + 1) * I) \ B

    if isnothing(m)
        m = select_nof_abscissas_gl(κ^(1 / 4), ϵ / 2)
    end

    G1, G2 = zero(A * B), zero(A * B)
    u, w = gausslegendre(m)
    for k in 1:m
        G1 += w[k] * F1(u[k])
        G2 += w[k] * F2(u[k])
    end
    G = ((cc - 1) * c * A - I) * G1 - (-c * A + (cc - 1) * I) * G2 - log(c) * B
    return LogmGLResult(G, 2m)
end

struct LogmDEResult
    X
    m
    h
    l
    r
end

"""
Computing log(A)B using the double exponential formula.
"""
function logm_de(A, λ_max, λ_min; B=I, ϵ=1e-12, m::Union{Nothing,Int}=nothing)
    κ = λ_max / λ_min
    c = 1 / sqrt(λ_max * λ_min)
    F(t) = ((1 + t) * c * A + (1 - t) * I) \ B

    l, r = get_interval(κ^(1 / 2), κ^(-1 / 2), ϵ / 2)
    if isnothing(m)
        m = select_nof_abscissas_de(κ^(1 / 2), l, r, ϵ / 2)
    end
    x = LinRange(l, r, m)
    h = (r - l) / (m - 1)
    w = fill(h, m)
    w[1], w[m] = h / 2, h / 2
    t = @. tanh(π / 2 * sinh(x))
    dt = @. π / 2 * cosh(x) * sech(π / 2 * sinh(x))^2

    Tl, Tr = zero(A * B), zero(A * B)
    for k in 1:(m ÷ 2)
        Tl += w[k] * dt[k] * F(t[k])
    end
    for k in m:-1:(m ÷ 2 + 1)
        Tr += w[k] * dt[k] * F(t[k])
    end
    T = (c * A - I) * (Tl + Tr) - log(c) * B
    return LogmDEResult(T, m, h, l, r)
end

function get_interval(λ_max, λ_min, ϵ)
    norm_A_minus_I = λ_max - 1
    norm_A_inv = 1 / λ_min
    ϵ_max = 3 * norm_A_minus_I * norm_A_inv / (1 + norm_A_inv)
    ϵ = ϵ < ϵ_max ? ϵ : ϵ_max / 2

    a = min(ϵ / 3 / norm_A_minus_I, 1 / 2 / norm_A_minus_I)
    α = (log(a) - log1p(-a)) / 2 # = atanh(2a-1)
    l = asinh(2 * α / pi)
    
    δ = a / norm_A_inv
    b1 = 1 - δ
    b2 = 2 * norm_A_inv / (2 * norm_A_inv + 1)
    if b1 >= b2
        β = (log1p(-δ) - log(δ))/2 # = atanh(2b1-1)
    else
        β = (log(b2) - log1p(-b2)) / 2 # = atanh(2b2-1)
    end
    r = asinh(2 * β / pi)
    return l, r
end

function select_nof_abscissas_de(λ, l, r, ϵ)
    function f(x)
        return (λ - 1) * π / 2 * cosh(x) * sech(π / 2 * sinh(x))^2 /
               ((1 + tanh(π / 2 * sinh(x))) * λ + (1 - tanh(π / 2 * sinh(x))))
    end
    err = Inf
    exact = log(λ)
    m_max = 1000
    m_min = 2
    while m_max - m_min > 1
        m = (m_max + m_min) ÷ 2
        x = LinRange(l, r, m)
        h = (r - l) / (m - 1)
        w = fill(h, m)
        w[1], w[m] = h / 2, h / 2
        approx = sum(@. w * f(x))
        err = abs(exact - approx)
        if err > ϵ
            m_min = m
        else
            m_max = m
        end
    end
    return m_max
end

"""
Computing log(A)B using the preconditioned double exponential formula.
"""
function logm_pde(A, λ_max, λ_min; B=I, ϵ=1e-12, m::Union{Nothing,Int}=nothing)
    κ = λ_max / λ_min
    c = 1 / sqrt(λ_max * λ_min)
    cc = sqrt((sqrt(κ) + 1) * (1 / sqrt(κ) + 1)) # = c'
    F1(t) = (((cc - 1) * t + 1 + cc) * c * A + (1 - t) * I) \ B
    F2(t) = ((1 - t) * c * A + ((cc - 1) * t + cc + 1) * I) \ B

    l, r = get_interval(κ^(1 / 4), κ^(-1 / 4), ϵ / 2)
    m = select_nof_abscissas_de(κ^(1 / 4), l, r, ϵ / 2)
    x = LinRange(l, r, m)
    h = (r - l) / (m - 1)
    w = fill(h, m)
    w[1], w[m] = h / 2, h / 2
    t = @. tanh(π / 2 * sinh(x))
    dt = @. π / 2 * cosh(x) * sech(π / 2 * sinh(x))^2

    T1l, T1r, T2l, T2r = zero(A * B), zero(A * B), zero(A * B), zero(A * B)
    for k in 1:(m ÷ 2)
        T1l += w[k] * dt[k] * F1(t[k])
        T2l += w[k] * dt[k] * F2(t[k])
    end
    for k in m:-1:(m ÷ 2 + 1)
        T1r += w[k] * dt[k] * F1(t[k])
        T2r += w[k] * dt[k] * F2(t[k])
    end
    T =
        ((cc - 1) * c * A - I) * (T1l + T1r) - (-c * A + (cc - 1) * I) * (T2l + T2r) -
        log(c) * B
    return LogmDEResult(T, 2m, h, l, r)
end

end
