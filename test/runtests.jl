using SymTridiagonalQL
using Test

using LinearAlgebra
using Random
import SymTridiagonalQL: eigvals_2x2, eigvals_3x3, dotu

eigperm(λ) = sortperm(λ, by=z->(abs(z),angle(z)))

function compare_eigvals(λ, λref)
    # This is not a fool-proof way of sorting the eigenvalues for
    # comparison, but good enough for our test cases.
    p1 = eigperm(λref)
    p2 = eigperm(λ)

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

function test_eigenvectors(H, ee)
    # Test if the purported eigenvectors of H indeed are such
    for (i,λ) in enumerate(ee.values)
        ϕ = view(ee.vectors, :, i)
        H*ϕ ≈ λ*ϕ || @error "Eigenpair #$(i) with λ = $(λ) is, in fact, not an eigenpair"
        @test H*ϕ ≈ λ*ϕ
    end
end

function compare_eigen(ee, eeref)
    λref = eeref.values
    λ = ee.values

    p1 = eigperm(λref)
    p2 = eigperm(λ)

    @test λref[p1] ≈ λ[p2]

    for (i,j) in zip(p1,p2)
        ϕref = view(eeref.vectors, :, i)
        ϕ = view(ee.vectors, :, j)

        abs(dotu(ϕref,ϕ)) ≈ 1 ||
            @error "Eigenvector with λ = $(λ[j]) does not agree with reference eigenvector with λ = $(λref[i])"

        @test abs(dotu(ϕref,ϕ)) ≈ 1
    end

    # Also test orthonormality by rotating each column of the overlap
    # matrix, such that the diagonal elements are purely real,
    # positive.
    S = complex(transpose(eeref.vectors[:,p1])*ee.vectors[:,p2])
    for j = 1:size(S,1)
        s = S[j,j]
        lmul!(exp(-im*angle(s)), view(S, :, j))
    end
    if size(S) == (200,200)
        @test S ≈ I rtol=1e-6
    else
        @test S ≈ I
    end
end

function compare_eigen(H::SymTridiagonal)
    ee = triql_vectors(H)
    test_eigenvectors(H, ee)

    eeref = eltype(H) <: Real ? eigen(real(H)) : eigen(Matrix(H))
    Q = eeref.vectors
    for j in 1:size(Q,2)
        q = view(Q, :, j)
        n = dotu(q,q)
        q ./= √n
    end

    compare_eigen(ee, eeref)
end

@testset "SymTridiagonalQL.jl" begin
    Random.seed!(1234)

    @testset "$(n)×$(n) matrix eigenvalues/eigenpairs" for n = [1, 2, 3, 4,10,100]
        for (D,E,degenerate) = [(-2ones(n), ones(n-1), false),
                                (-2ones(2n), vcat(ones(n-1),0,ones(n-1)), true),
                                (complex(-2ones(n)), complex(ones(n-1)), false),
                                (collect(-2*(1:n) .+ 0.4im), collect(0.4*(1:n-1) .+ 0im), false),
                                (collect(1.0:2n), vcat(1.0:(n-1),0,1:(n-1)), false),
                                (vcat(1.0:n, 1.0:n), vcat(1.0:(n-1),0,1:(n-1)), true),
                                (ones(n), zeros(n-1), true),
                                (collect(1.0:n), zeros(n-1), false),
                                (collect(1.0:n), ones(n-1) .+ 3, false),
                                (collect(1.0:n), n > 1 ? vcat(ones(n-2),0.0) : Float64[], false),
                                (collect(1.0:n), n > 1 ? vcat(0.0,ones(n-2)) : Float64[], false),
                                (collect(1.0:n)*im, zeros(n-1)*im, false),
                                (collect(1.0:n)*im, 0.1 .+ zeros(n-1)*im, false),
                                (vcat(ones(n-1),0), zeros(n-1), true),
                                (rand(n), rand(n-1), false),
                                (rand(ComplexF64, n), rand(ComplexF64, n-1), false)]
            H = SymTridiagonal(D, E)

            compare_eigvals(H)

            # We cannot compare the eigenvectors of degenerate spectra
            # easily, since eigenvectors belonging to the same
            # eigenvalue, may appear in any order.
            degenerate || compare_eigen(H)
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

    include("rotations.jl")
end
