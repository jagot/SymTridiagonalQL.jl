function wilkinson_shift2(Dᵢ, Dᵢ₊₁, Eᵢ; verbosity=0)
    iszero(Eᵢ) && return zero(Dᵢ)

    p = (Dᵢ₊₁ - Dᵢ)/2Eᵢ
    q = √(p^2 + 1)

    a = p+q
    b = p-q
    c = abs(a) < abs(b) ? a : b

    σ = Dᵢ + Eᵢ*c

    if verbosity > 0
        σ₊ = Dᵢ + Eᵢ*a
        σ₋ = Dᵢ + Eᵢ*b
        @info "Wilkinson 2×2 shift" Dᵢ Dᵢ₊₁ Eᵢ σ₊ σ₋ σ
    end

    σ
end

function wilkinson_shift2(D::AbstractVector{T}, E::AbstractVector, i; kwargs...) where T
    n = length(D)
    n == 0 && return zero(T)
    n == 1 && return D[1]
    wilkinson_shift2(D[i], D[i+1], E[i]; kwargs...)
end

wilkinson_shift2(T::SymTridiagonal, i; kwargs...) =
    wilkinson_shift2(T.dv, T.ev; kwargs...)

function wilkinson_shift3(Dᵢ, Dᵢ₊₁, Dᵢ₊₂, Eᵢ, Eᵢ₊₁; verbosity=0)
    iszero(Eᵢ₊₁) && return wilkinson_shift2(Dᵢ, Dᵢ₊₁, Eᵢ; verbosity=verbosity)
    a = (Dᵢ + Dᵢ₊₁ + Dᵢ₊₂)

    Y = -a^2 +
        3*(Dᵢ₊₁*Dᵢ₊₂ + Dᵢ*(Dᵢ₊₁ + Dᵢ₊₂) - Eᵢ^2 - Eᵢ₊₁^2)

    b = (2Dᵢ - Dᵢ₊₁ - Dᵢ₊₂)
    Z = -(Dᵢ + Dᵢ₊₁ - 2Dᵢ₊₂)*(a*(Dᵢ - 2Dᵢ₊₁ + Dᵢ₊₂) + 9Eᵢ^2) +
        9a*Eᵢ₊₁^2

    a3 = a/3
    c = complex(Z + √(Z^2 + 4Y^3))^(1/3)
    d = 1 + im*√3
    e = ∛2

    Λ₁ = a3 + e*Y/3c - c/3e
    Λ₂ = a3 - d*Y/(3e*c) + conj(d)*c/6e
    Λ₃ = a3 - conj(d)*Y/(3e*c) + d*c/6e

    σ = Λ₁
    abs(Dᵢ - σ) > abs(Dᵢ - Λ₂) && (σ = Λ₂)
    abs(Dᵢ - σ) > abs(Dᵢ - Λ₃) && (σ = Λ₃)

    verbosity > 0 &&
        @info "Wilkinson 3×3 shift" Dᵢ Dᵢ₊₁ Dᵢ₊₂ Eᵢ Eᵢ₊₁ Y Z Λ₁ Λ₂ Λ₃ σ

    σ
end

function wilkinson_shift3(D::AbstractVector{T}, E::AbstractVector, i; kwargs...) where T
    n = length(D)
    n == 0 && return zero(T)
    n == 1 && return D[1]
    n == 2 && return wilkinson_shift2(D, E, i; kwargs...)
    wilkinson_shift3(D[i], D[i+1], D[i+2],
                     E[i], E[i+1]; kwargs...)
end

wilkinson_shift3(T::SymTridiagonal, i; kwargs...) =
    wilkinson_shift3(T.dv[i], T.dv[i+1], T.dv[i+2],
                     T.ev[i], T.ev[i+2]; kwargs...)
