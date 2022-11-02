using SymTridiagonalQL
using Test

using LinearAlgebra

@testset "SymTridiagonalQL.jl" begin
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
