using SymTridiagonalQL
using Test

using LinearAlgebra
using Random
import SymTridiagonalQL: eigvals_2x2, eigvals_3x3

function compare_eigvals(λ, λref)
    # This is not a fool-proof way of sorting the eigenvalues for
    # comparison, but good enough for our test cases.
    p1 = sortperm(λref, by=z->(abs(z),angle(z)))
    p2 = sortperm(λ, by=z->(abs(z),angle(z)))

    @test λref[p1] ≈ λ[p2]
end

function compare_eigvals(H::SymTridiagonal)
    λ = diag(triql(H))
    λref = isreal(H) ? eigvals(real(H)) : eigvals(Matrix(H))
    compare_eigvals(λ, λref)

    # Additionally, for small matrices, we also compare the explicit
    # formulas we have implemented for the shifts:
    D = H.dv
    E = H.ev
    n = length(D)
    if n == 2
        compare_eigvals(collect(eigvals_2x2(D[1], D[2], E[1])), λref)
    elseif n == 3
        compare_eigvals(collect(eigvals_3x3(D[1], D[2], D[3], E[1], E[2], verbosity=Inf)), λref)
    end
end

@testset "SymTridiagonalQL.jl" begin
    Random.seed!(1234)

    @testset "$(n)×$(n) matrix eigenvalues" for n = [1,2,3,4,10,100]
        for (D,E) = [(-2ones(n), ones(n-1)),
                     (-2ones(2n), vcat(ones(n-1),0,ones(n-1))),
                     (complex(-2ones(n)), complex(ones(n-1))),
                     (collect(-2*(1:n) .+ 0.4im), collect(0.4*(1:n-1) .+ 0im)),
                     (collect(1.0:2n), vcat(1.0:(n-1),0,1:(n-1))),
                     (vcat(1.0:n, 1.0:n), vcat(1.0:(n-1),0,1:(n-1))),
                     (ones(n), zeros(n-1)),
                     (collect(1.0:n), zeros(n-1)),
                     (collect(1.0:n), ones(n-1) .+ 3),
                     (collect(1.0:n), n > 1 ? vcat(ones(n-2),0.0) : Float64[]),
                     (collect(1.0:n), n > 1 ? vcat(0.0,ones(n-2)) : Float64[]),
                     (collect(1.0:n)*im, zeros(n-1)*im),
                     (collect(1.0:n)*im, 0.1 .+ zeros(n-1)*im),
                     (vcat(ones(n-1),0), zeros(n-1)),
                     (rand(n), rand(n-1)),
                     (rand(ComplexF64, n), rand(ComplexF64, n-1))]
            H = SymTridiagonal(D, E)

            compare_eigvals(H)
        end
    end

    @testset "Diagonal matrix, n = $n" for n = [1,2,3,4,10,40]
        λref = collect(1.0:n)
        H = SymTridiagonal(λref, zeros(n-1))
        @testset "Shift-mode = $(shift_mode), verbosity = $(verbosity)" for
            (shift_mode,verbosity) in [(3,0),(3,4),(2,0),(2,4)]
            λ = diag(triql(H, verbosity=verbosity, shift_mode=shift_mode))
            @test λ ≈ λref
        end
    end
end
