#=

Equations below refer to

- Noble, J., Lubasch, M., Stevens, J., & Jentschura,
  U. (2017). Diagonalization of complex symmetric matrices:
  generalized Householder reflections, iterative deflation and
  implicit shifts. Computer Physics Communications, 221(nil),
  304–316. http://dx.doi.org/10.1016/j.cpc.2017.06.014

=#

function eigvals_2x2(Dᵢ, Dᵢ₊₁, Eᵢ)
    iszero(Eᵢ) && return Dᵢ, Dᵢ₊₁

    # Eq. (14)

    p = (Dᵢ₊₁ - Dᵢ)/2Eᵢ
    q = √(p^2 + 1)

    a = p+q
    b = p-q

    Λ₊ = Dᵢ + Eᵢ*a
    Λ₋ = Dᵢ + Eᵢ*b

    Λ₊, Λ₋
end

function wilkinson_shift(Dᵢ, Dᵢ₊₁, Eᵢ; verbosity=0)
    iszero(Eᵢ) && return zero(Dᵢ)

    # We could early-out on computing only the necessary eigenvalue
    # using
    #
    # c = abs(a) < abs(b) ? a : b
    #
    # but we opt for clarity. This is hardly a bottleneck, but could
    # be optimized down the road by passing some kind of extra
    # argument to eigvals_2x2.

    σ₊,σ₋ = eigvals_2x2(Dᵢ, Dᵢ₊₁, Eᵢ)

    σ = abs(σ₊-Dᵢ) < abs(σ₋-Dᵢ) ? σ₊ : σ₋

    verbosity > 0 &&
        @info "Wilkinson 2×2 shift" Dᵢ Dᵢ₊₁ Eᵢ σ₊ σ₋ σ

    σ
end

function wilkinson_shift(D::AbstractVector{T}, E::AbstractVector, i; kwargs...) where T
    n = length(D)
    n == 0 && return zero(T)
    n == 1 && return D[1]
    wilkinson_shift(D[i], D[i+1], E[i]; kwargs...)
end

wilkinson_shift(T::SymTridiagonal, i; kwargs...) =
    wilkinson_shift(T.dv, T.ev, i; kwargs...)

function eigvals_3x3(Dᵢ, Dᵢ₊₁, Dᵢ₊₂, Eᵢ, Eᵢ₊₁; verbosity=0)
    if iszero(Eᵢ)
        return (Dᵢ, eigvals_2x2(Dᵢ₊₁, Dᵢ₊₂, Eᵢ₊₁)...)
    elseif iszero(Eᵢ₊₁)
        return (Dᵢ₊₂, eigvals_2x2(Dᵢ, Dᵢ₊₁, Eᵢ)...)
    end

    # Eq. (16); there appears to be a typo in Eqs. (16b) & (16c),
    # where 2^(1/3) should be replaced by 2^(2/3) in the second term
    # of each of the quoted equations. This was figured out using by
    # asking WolframAlpha for "eigenvalues of
    # {{a,b,0},{b,d,e},{0,e,f}}":
    # https://www.wolframalpha.com/input?i=eigenvalues+of+%7B%7Ba%2Cb%2C0%7D%2C%7Bb%2Cd%2Ce%7D%2C%7B0%2Ce%2Cf%7D%7D

    a = (Dᵢ + Dᵢ₊₁ + Dᵢ₊₂)
    a3 = a/3

    Y = -a^2 +
        3*(Dᵢ₊₁*Dᵢ₊₂ + Dᵢ*(Dᵢ₊₁ + Dᵢ₊₂) - Eᵢ^2 - Eᵢ₊₁^2)

    b = (2Dᵢ - Dᵢ₊₁ - Dᵢ₊₂)
    Z = -(Dᵢ + Dᵢ₊₁ - 2Dᵢ₊₂)*(b*(Dᵢ - 2Dᵢ₊₁ + Dᵢ₊₂) + 9Eᵢ^2) +
        9b*Eᵢ₊₁^2

    c = (Z + √(complex(Z^2 + 4Y^3)))^(1/3)
    d = 1 + im*√3
    e = ∛2

    A₁ = e*Y/3c
    B₁ = -c/3e

    A₂ = -d*Y/(3e^2*c)
    B₂ = conj(d)*c/6e

    A₃ = -conj(d)*Y/(3e^2*c)
    B₃ = d*c/6e

    Λ₁ = a3 + A₁ + B₁
    Λ₂ = a3 + A₂ + B₂
    Λ₃ = a3 + A₃ + B₃

    verbosity > 0 &&
        @info "3×3 eigenvalues" Dᵢ Dᵢ₊₁ Dᵢ₊₂ Eᵢ Eᵢ₊₁ a a3 b Y Z c d e A₁ B₁ A₂ B₂ A₃ B₃ Λ₁ Λ₂ Λ₃

    Λ₁,Λ₂,Λ₃
end

function cubic_shift(Dᵢ, Dᵢ₊₁, Dᵢ₊₂, Eᵢ, Eᵢ₊₁; verbosity=0)
    iszero(Eᵢ₊₁) && return wilkinson_shift(Dᵢ, Dᵢ₊₁, Eᵢ; verbosity=verbosity)

    Λ₁,Λ₂,Λ₃ = eigvals_3x3(Dᵢ, Dᵢ₊₁, Dᵢ₊₂, Eᵢ, Eᵢ₊₁; verbosity=verbosity)

    σ = Λ₁
    abs(Dᵢ - σ) > abs(Dᵢ - Λ₂) && (σ = Λ₂)
    abs(Dᵢ - σ) > abs(Dᵢ - Λ₃) && (σ = Λ₃)

    verbosity > 0 &&
        @info "Cubic 3×3 shift" Λ₁ Λ₂ Λ₃ σ

    σ
end

function cubic_shift(D::AbstractVector{T}, E::AbstractVector, i; kwargs...) where T
    length(D) < 3 && return wilkinson_shift(D, E, i; kwargs...)

    cubic_shift(D[i], D[i+1], D[i+2],
                E[i], E[i+1]; kwargs...)
end

cubic_shift(T::SymTridiagonal, i; kwargs...) =
    cubic_shift(T.dv, T.ev, i; kwargs...)

function get_shift(D, E, shift_mode; kwargs...)
    if shift_mode == 0
        zero(eltype(D))
    elseif shift_mode == 1
        D[1]
    elseif shift_mode == 2
        wilkinson_shift(D, E, 1; kwargs...)
    elseif shift_mode == 3
        cubic_shift(D, E, 1; kwargs...)
    else
        throw(ArgumentError("Invalid shift_mode, choose 0–3"))
    end
end
