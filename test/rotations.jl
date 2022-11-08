import SymTridiagonalQL: Rotation, CompositeRotation

@testset "Rotations" begin
    R = Rotation(π/3, 2)
    R2 = Rotation(π/4, 3)

    @test_throws DimensionMismatch Matrix(R, (1,1))

    M = Matrix(R, (5,5))
    @test M[2,2] == -R.c
    @test M[2,3] == R.s
    @test M[3,2] == R.s
    @test M[3,3] == R.c
    @test M^2 ≈ I

    @test lmul!(R, Matrix(1.0I, (5,5))) ≈ M

    @test lmul!(R*R2, Matrix(1.0I, (5,5))) ≈ M*Matrix(R2, (5,5))
end
