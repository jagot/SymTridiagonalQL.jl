macro check_E(E,i,ϵ,verbosity)
    quote
        let i = $(esc(i))
            if i > 1 && abs($(esc(E))[i]) < $(esc(ϵ))
                $(esc(verbosity)) > 0 && @info "Premature zero at E[$i]"
                return i
            end
        end
    end
end

function chase_the_bulge!(D, E, c, s; verbosity=0)
    n = length(D)
    n == 1 && return

    fmt = nothing
    fmt′ = nothing
    if verbosity > 0
        nd = length(digits(n))
        @info "Chasing the bulge" n c s c^2 + s^2
        printfmtln("{1:>$(nd)s}  | {2:22s} {3:22s} | {4:22s} {5:22s} {6:22s} {7:22s} {8:22s} | {9:22s}",
                   "i", "c", "s",
                   "Dᵢ₊₁", "Dᵢ", "Eᵢ₊₁", "Eᵢ", "Eᵢ₋₁", "F")

        fmt = FormatExpr("{1:>$(nd)d}  | {2:+0.3e}{3:+0.3e}im {4:+0.3e}{5:+0.3e}im | {6:+0.3e}{7:+0.3e}im {8:+0.3e}{9:+0.3e}im {10:+0.3e}{11:+0.3e}im {12:+0.3e}{13:+0.3e}im {14:+0.3e}{15:+0.3e}im | {16:+0.3e}{17:+0.3e}im")
        fmt′ = FormatExpr("{1:>$(nd)d}′ | {2:+0.3e}{3:+0.3e}im {4:+0.3e}{5:+0.3e}im | {6:+0.3e}{7:+0.3e}im {8:+0.3e}{9:+0.3e}im {10:+0.3e}{11:+0.3e}im {12:+0.3e}{13:+0.3e}im {14:+0.3e}{15:+0.3e}im | {16:+0.3e}{17:+0.3e}im")
    end

    Dᵢ₊₁ = D[end]
    Dᵢ = D[end-1]
    Eᵢ = E[end]
    Eᵢ₋₁ = n > 2 ? E[end-1] : zero(Eᵢ)
    Eᵢ₊₁ = zero(Eᵢ)
    F = zero(Eᵢ)

    Dᵢ₊₁′ = zero(Dᵢ₊₁)
    Dᵢ′ = zero(Dᵢ)
    Eᵢ₋₁′ = zero(Eᵢ₋₁)
    Eᵢ′ = zero(Eᵢ)
    Eᵢ₊₁′ = zero(Eᵢ₊₁)
    F′ = zero(F)

    c′ = zero(c)
    s′ = zero(s)

    ϵ = eps(real(eltype(E)))
    # We estimate the norm as ||A||₂ ≤ 3 maxᵢ,ⱼ |aᵢⱼ| since ||A||₂ ≤
    # √(||A||₁||A||_∞) by Schur's test, see
    # https://math.stackexchange.com/a/3743734/45104
    N = 3max(maximum(abs, D),maximum(abs, E))
    tol = √ϵ*N

    verbosity > 0 && @show ϵ N tol norm(SymTridiagonal(D, E))

    for i = n-1:-1:1
        if i < n-1
            B = inv(√(Eᵢ₊₁^2 + F^2))
            c = Eᵢ₊₁*B
            s = F*B
        end

        c² = c*c
        cs = c*s
        s² = s*s

        a = 2cs*Eᵢ

        Dᵢ₊₁′ = D[i+1] = c²*Dᵢ₊₁ + a + s²*Dᵢ
        Dᵢ′ = D[i] = c²*Dᵢ - a + s²*Dᵢ₊₁
        if i < n-1
            Eᵢ₊₁′ = E[i+1] = √(Eᵢ₊₁^2+F^2)
            # @check_E E i+1 ϵ verbosity
        end
        Eᵢ′ = E[i] = (c² - s²)*Eᵢ + cs*(Dᵢ - Dᵢ₊₁)
        @check_E E i ϵ verbosity
        if i > 1
            Eᵢ₋₁′ = E[i-1] = c*Eᵢ₋₁
            @check_E E i-1 ϵ verbosity
            F′ = s*Eᵢ₋₁
        end

        if verbosity > 0
            printfmtln(fmt, i,
                       real(c), imag(c),
                       real(s), imag(s),
                       real(Dᵢ₊₁), imag(Dᵢ₊₁),
                       real(Dᵢ), imag(Dᵢ),
                       real(Eᵢ₊₁), imag(Eᵢ₊₁),
                       real(Eᵢ), imag(Eᵢ),
                       real(Eᵢ₋₁), imag(Eᵢ₋₁),
                       real(F), imag(F))
            printfmtln(fmt′, i,
                       real(c′), imag(c′),
                       real(s′), imag(s′),
                       real(Dᵢ₊₁′), imag(Dᵢ₊₁′),
                       real(Dᵢ′), imag(Dᵢ′),
                       real(Eᵢ₊₁′), imag(Eᵢ₊₁′),
                       real(Eᵢ′), imag(Eᵢ′),
                       real(Eᵢ₋₁′), imag(Eᵢ₋₁′),
                       real(F′), imag(F′))
        end

        if i > 1
            Dᵢ₊₁ = Dᵢ′
            Dᵢ = D[i-1]
            Eᵢ₊₁ = Eᵢ′
            Eᵢ = Eᵢ₋₁′
            Eᵢ₋₁ = E[i-1]
            F = F′
            # c = c′
            # s = s′
        end
    end
    1
end

function triql!(D, E; max_iter=150, tol=100eps(real(eltype(E))), verbosity=0, shift_mode=2)
    nD = length(D)
    nD < 2 && return 1
    nD < 3 && (shift_mode = min(shift_mode, 2))

    fmt = if verbosity > 0
        nd = length(digits(max_iter))
        printfmtln("{1:>$(nd)s} {2:$(nd)s} | {3:10s} {4:10s}",
                   "n", "", "|Dⱼ|", "|Eⱼ|")
        FormatExpr("{1:>$(nd)d}/{2:$(nd)d} | {3:+0.3e} {4:+0.3e}")
    end


    for n = 1:max_iter
        σ = get_shift(D, E, shift_mode; verbosity=verbosity-2)

        Dₙ = D[end]
        Eₙ₋₁ = E[end]
        iszero(Eₙ₋₁) && return length(E)
        A = inv(√((Dₙ - σ)^2 + Eₙ₋₁^2))
        c = (Dₙ - σ)*A
        s = Eₙ₋₁*A
        verbosity > 1 && @info "Initial rotation" Dₙ Eₙ₋₁ A c s c^2 + s^2

        i = chase_the_bulge!(D, E, c, s, verbosity=verbosity-2)
        if verbosity > 0
            verbosity > 1 && display(SymTridiagonal(D, E))
            printfmtln(fmt, n, max_iter, abs(D[1]), abs(E[1]))
        end

        i ≠ 1 && return i
        abs(E[1]) < tol && return 1
    end

    abs(E[1]) < tol || @warn "Convergence not reached in $(max_iter) iterations" tol
    1
end

function triql!(T::SymTridiagonal; verbosity=0, kwargs...)
    stack = [1:size(T,1)]
    while !isempty(stack)
        r = pop!(stack)
        n = r[end]
        verbosity > 0 && @show r
        for j ∈ r
            verbosity > 0 && @show j
            i = triql!(view(T.dv, j:n), view(T.ev, j:n-1);
                       verbosity=verbosity, kwargs...)
            if i ≠ 1
                ar = r[1]:r[i]
                br = r[i+1]:n
                push!(stack, ar)
                push!(stack, br)

                verbosity > 0 &&
                    @info "Premature zero at j = $(r[i]), deflating into $(ar) and $(br)"
                break
            end
            # @warn "Breaking"
            # break
        end
    end
    T
end

triql(T; kwargs...) = triql!(copy(T); kwargs...)

export triql!, triql
