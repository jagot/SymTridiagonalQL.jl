import LinearAlgebra: sorteig!

#=

References:

- Parlett, B. N. (1998). The Symmetric Eigenvalue Problem. : Society
  for Industrial and Applied
  Mathematics. https://doi.org/10.1137/1.9781611971163

- Cullum, J. K., & Willoughby, R. A. (1996). A $QL$ procedure for
  computing the eigenvalues of complex symmetric tridiagonal
  matrices. SIAM Journal on Matrix Analysis and Applications, 17(1),
  83–109. http://dx.doi.org/10.1137/s0895479894137639

=#

simple_test(Eᵢ, tol) = abs(Eᵢ) < tol

function expensive_test(Eᵢ, Eᵢ₋₁, Eᵢ₊₁, Dᵢ, Dᵢ₊₁, ϵ)
    # Parlett (1998), exercise 7.11.4
    δ = (Dᵢ₊₁ - Dᵢ)/2 # Eq. (7.24)
    # We need to expliclity handle the case when two successive
    # diagonal elements are equal and the off-diagonal element is
    # zero.
    iszero(δ) && iszero(Eᵢ) ||
        abs(Eᵢ)*(abs(Eᵢ₋₁) + abs(Eᵢ₊₁)) < ϵ*abs(δ)*(abs(Dᵢ) + abs(Dᵢ₊₁))
end

is_small(Eᵢ, Eᵢ₋₁, Eᵢ₊₁, Dᵢ, Dᵢ₊₁, tol, ϵ) =
    iszero(Eᵢ) ||
    simple_test(Eᵢ, tol) &&
    expensive_test(Eᵢ, Eᵢ₋₁, Eᵢ₊₁, Dᵢ, Dᵢ₊₁, ϵ)

function findnext_small(E, D, l, tol, ϵ; verbose=false)
    N = length(D)
    Dᵢ = D[l]
    Eᵢ = E[l]
    Eᵢ₋₁ = l > 1 ? E[l-1] : zero(Eᵢ)
    for i = l:N
        Dᵢ₊₁ = i < N ? D[i+1] : zero(Dᵢ)
        Eᵢ₊₁ = i < N ? E[i+1] : zero(Eᵢ)
        verbose &&
            @info "" i (Eᵢ, Eᵢ₋₁, Eᵢ₊₁) (Dᵢ, Dᵢ₊₁) simple_test(Eᵢ, tol) expensive_test(Eᵢ, Eᵢ₋₁, Eᵢ₊₁, Dᵢ, Dᵢ₊₁, ϵ) is_small(Eᵢ, Eᵢ₋₁, Eᵢ₊₁, Dᵢ, Dᵢ₊₁, tol, ϵ)
        is_small(Eᵢ, Eᵢ₋₁, Eᵢ₊₁, Dᵢ, Dᵢ₊₁, tol, ϵ) && return i

        Dᵢ = Dᵢ₊₁

        Eᵢ₋₁ = Eᵢ
        Eᵢ = Eᵢ₊₁
    end
    nothing
end

# We estimate the norm as ||A||₂ ≤ 3 maxᵢ,ⱼ |aᵢⱼ| since ||A||₂ ≤
# √(||A||₁||A||_∞) by Schur's test, see
# https://math.stackexchange.com/a/3743734/45104
estimate_tridiag_norm(D, E) = 3max(maximum(abs, D),maximum(abs, E))

macro check_E(D,E,i,tol,ϵ,verbosity)
    quote
        let i = $(esc(i)), E = $(esc(E)), D = $(esc(D))
            Eᵢ = E[i]
            # We follow section 7.11.3 of Parlett (1998); first we
            # perform a simple test, to see if the off-diagonal
            # element magnitude |Eᵢ| is larger than a set
            # tolerance. If not, we perform a more expensive test, to
            # see if we should consider the element to be zero.
            if i > 1 && simple_test(Eᵢ, $(esc(tol)))
                $(esc(verbosity)) > 0 && @info "Possibly premature zero at E[$i], investigating"

                n = length(D)
                Eᵢ₊₁ = i == n-1 ? zero(Eᵢ) : E[i+1]
                if expensive_test(Eᵢ, E[i-1], Eᵢ₊₁, D[i], D[i+1], $(esc(ϵ)))
                    $(esc(verbosity)) > 0 && @info "Premature zero at E[$i]"
                    E[i] = 0
                    return i
                end
            end
        end
    end
end

function cancellation_imminent(w, ι=100eps(real(w)))
    a = real(w)
    b = imag(w)
    1 + a^2 ≤ (1+ι)*b^2
end

@doc raw"""
    cullum1996_sweep!(D, E, shift_mode, verbosity)

Implementation of the CMTQL1 algorithm, as described in section 5.3 of

- Cullum, J. K., & Willoughby, R. A. (1996). A $QL$ procedure for
  computing the eigenvalues of complex symmetric tridiagonal
  matrices. SIAM Journal on Matrix Analysis and Applications, 17(1),
  83–109. http://dx.doi.org/10.1137/s0895479894137639

Specifically, this covers steps 3–4, since the other steps are taken
care of by the driver routine [`triql!`](@ref).

"""
function cullum1996_sweep!(k, D, E, shift_mode, verbosity; Qright=nothing, max_sweep_tries=10)
    U = promote_type(eltype(D), eltype(E))

    M = length(D)

    for j = 1:max_sweep_tries
        ϕ = (k == 0 ?
            rand(promote_type(eltype(D), eltype(E)))*estimate_tridiag_norm(D, E) :
            get_shift(D, E, shift_mode, verbosity=verbosity-2))

        s = one(U)
        c = -one(U)
        w = one(U)
        F = zero(U)
        G = D[end] - ϕ
        B = zero(U)
        ϕ = zero(U)

        fmt = nothing
        fmt′ = nothing
        if verbosity > 0
            @info "Cullum & Willoughby sweep" j k ϕ M c s c^2 + s^2
            nd = length(digits(M))
            printfmtln("{1:>$(nd)s} | {2:22s} {3:22s} | {4:22s} {5:22s} {6:22s} {7:22s} | {8:22s} {9:22s}",
                       "i", "c", "s",
                       "Dᵢ₊₁", "Dᵢ", "Eᵢ₊₁", "Eᵢ", "F", "G")

            fmt = FormatExpr("{1:>$(nd)d} | {2:+0.3e}{3:+0.3e}im {4:+0.3e}{5:+0.3e}im | {6:+0.3e}{7:+0.3e}im {8:+0.3e}{9:+0.3e}im {10:+0.3e}{11:+0.3e}im {12:+0.3e}{13:+0.3e}im | {14:+0.3e}{15:+0.3e}im {16:+0.3e}{17:+0.3e}im")
            fmt′ = FormatExpr("{1:>$(nd)s} | {2:+0.3e}{3:+0.3e}im {4:+0.3e}{5:+0.3e}im")
        end

        for i = M-1:-1:1
            # This is lagging one iteration behind, i.e. these are the
            # values F and E[i] would have from the previous rotation.
            F = s*E[i]
            B = -c*E[i]

            verbosity > 0 &&
                printfmtln(fmt, i,
                           real(c), imag(c),
                           real(s), imag(s),
                           real(D[i+1]), imag(D[i+1]),
                           real(D[i]), imag(D[i]),
                           real(E[i+1]), imag(E[i+1]),
                           real(E[i]), imag(E[i]),
                           real(F), imag(F),
                           real(G), imag(G))

            E[i+1] = if abs(G) > abs(F)
                w = F/G
                r = √(1 + w^2)
                c = inv(r)
                s = w*c
                r*G
            else
                w = G/F
                r = √(1 + w^2)
                s = inv(r)
                c = w*s
                r*F
            end
            verbosity > 0 &&
                printfmtln(fmt′, "",
                           real(c), imag(c),
                           real(s), imag(s))

            if cancellation_imminent(w)
                @warn "Cancellation may occur, restarting sweep" w
                k = 1
                break
            end

            verbosity > 1 && @show D[i], D[i+1]
            G = D[i+1] - ϕ
            rr = (D[i] - G)*s + 2c*B
            verbosity > 1 && @show ϕ, G, rr
            ϕ = s*rr
            verbosity > 1 && @show ϕ
            D[i+1] = G + ϕ
            verbosity > 1 && @show D[i], D[i+1]
            G = B - c*rr
            verbosity > 1 && @show G

            if !isnothing(Qright)
                if verbosity > 1
                    @info "Forming eigenvectors" i c s
                    @show Qright[:,i]
                    @show Qright[:,i+1]
                end

                qi = view(Qright, :, i)
                qip1 = view(Qright, :, i+1)
                for k = eachindex(qi)
                    u = qi[k]
                    v = qip1[k]
                    qi[k] = -c*u + s*v
                    qip1[k] = s*u + c*v
                    if verbosity > 1
                        @show k
                        display(parent(Qright))
                    end
                end
            end
        end
        D[1] -= ϕ
        E[1] = G
        E[M] = false
        return
    end
    @warn "Successful sweep not accomplished in $(max_sweep_tries) tries"
end

function triql!(D, E; Qright=nothing, max_iter=150, ϵ=eps(real(eltype(E))), verbosity=0, shift_mode=2, kwargs...)
    N = length(D)
    @assert length(E) == N
    N < 2 && return
    N < 3 && (shift_mode = min(shift_mode, 2))

    fmt = if verbosity > 0
        nd = length(digits(max_iter))
        printfmtln("{1:>$(nd)s} {2:$(nd)s} | {3:10s} {4:10s}",
                   "n", "", "|Dⱼ|", "|Eⱼ|")
        FormatExpr("{1:>$(nd)d}/{2:$(nd)d} | {3:+0.3e} {4:+0.3e}")
    end

    tol = √2*ϵ

    ϵ = eps(real(eltype(E)))
    small = √ϵ*estimate_tridiag_norm(D, E) # norm(T)

    for l = 1:N
        converged = false
        for n = 1:max_iter
            m = something(findnext_small(E, D, l, small, ϵ; verbose=verbosity>4), N)
            if m == l
                converged = true
                break
            end
            if verbosity > 0
                horizontal_line(color=:green)
                @info "Diagonalizing submatrix $(l)–$(m)"
            end

            sel = l:m
            Ds = view(D, sel)
            Es = view(E, sel)
            Qrights = isnothing(Qright) ? nothing : view(Qright, :, sel)

            cullum1996_sweep!(n, Ds, Es, shift_mode, verbosity-2; Qright=Qrights, kwargs...)

            if verbosity > 0
                verbosity > 1 && display(SymTridiagonal(D, E))
                printfmtln(fmt, n, max_iter, abs(Ds[1]), abs(Es[1]))
            end

            converged = abs(Es[1]) < tol*(abs(Ds[1])+abs(Ds[2]))
            converged && break
        end
        if verbosity > 5
            horizontal_line(color=:green)
            display(SymTridiagonal(D, E))
            horizontal_line(color=:green)
        end

        converged || @warn "Convergence not reached in $(max_iter) iterations" tol
    end
end

function triql!(T::SymTridiagonal; verbosity=0, kwargs...)
    N = size(T,1)
    D = T.dv
    E = T.ev

    resize!(E,N)
    E[N] = false

    triql!(D, E; verbosity=verbosity, kwargs...)

    resize!(E,N-1)

    T
end

triql(T; kwargs...) = triql!(copy(T); kwargs...)

function triql_vectors(T; kwargs...)
    U = eltype(T)
    Qright = Matrix(one(U)*I, size(T))
    T = triql!(copy(T); Qright=Qright, kwargs...)
    Eigen(sorteig!(diag(T), Qright)...)
end

export triql!, triql, triql_vectors
