struct Rotation{T}
    c::T
    s::T
    i::Int
end

function Rotation(θ, i::Int)
    s,c = sincos(θ)
    Rotation(c, s, i)
end

Base.show(io::IO, R::Rotation{T}) where T =
    write(io, "2×2 Rotation{$T} at ($(R.i),$(R.i))")

function Base.Matrix(R::Rotation{T}, (m,n)) where T
    i = R.i
    m > i && n > i ||
        throw(DimensionMismatch("Cannot instantiate $(R) in a matrix of size ($m,$n)"))

    M = Matrix(one(T)*I, (m,n))
    M[i,i] = -R.c
    M[i,i+1] = R.s
    M[i+1,i] = R.s
    M[i+1,i+1] = R.c
    M
end

function LinearAlgebra.lmul!(R::Rotation{T}, M::StridedArray{T}) where T
    c,s,i = R.c,R.s,R.i
    m = size(M,1)

    if 0 < i < m
        for Ipost in CartesianIndices(size(M)[2:end])
            a,b = M[i,Ipost],M[i+1,Ipost]
            M[i,Ipost] = -c*a + s*b
            M[i+1,Ipost] = s*a + c*b
        end
    elseif i == m
        for Ipost in CartesianIndices(size(M)[2:end])
            M[m,Ipost] *= -c
        end
    elseif i == 0
        for Ipost in CartesianIndices(size(M)[2:end])
            M[1,Ipost] *= c
        end
    end

    M
end

function Base.show(io::IO, mime::MIME"text/plain", R::Rotation)
    show(io, R)
    println(io, ":")
    M = [-R.c R.s; R.s R.c]
    Base.print_matrix(io, M)
end

struct CompositeRotation{T}
    # First element is applied first, when multiplying to the right
    rotations::Vector{Rotation{T}}
end

Base.show(io::IO, R::CompositeRotation{T}) where T =
    write(io, "$(length(R.rotations))-rotation CompositeRotation{$T}")

function Base.show(io::IO, mime::MIME"text/plain", R::CompositeRotation)
    show(io, R)
    isempty(R.rotations) || println(io, ":")
    for r in R.rotations
        show(io, mime, r)
        println(io)
    end
end

Base.one(::CompositeRotation{T}) where T=
    CompositeRotation(Vector{Rotation{T}}())

Base.one(::Type{CompositeRotation{T}}) where T=
    CompositeRotation(Vector{Rotation{T}}())

Base.copy(R::CompositeRotation) =
    CompositeRotation(copy(R.rotations))

function LinearAlgebra.rmul!(A::CompositeRotation{T}, B::Rotation{T}) where T
    pushfirst!(A.rotations, B)
    A
end

function LinearAlgebra.lmul!(A::Rotation{T}, B::CompositeRotation{T}) where T
    push!(B.rotations, A)
    B
end

function LinearAlgebra.rmul!(A::CompositeRotation{T}, B::CompositeRotation{T}) where T
    prepend!(A.rotations, B.rotations)
    A
end

function LinearAlgebra.lmul!(A::CompositeRotation{T}, B::CompositeRotation{T}) where T
    append!(B.rotations, A.rotations)
    B
end

Base.:(*)(A::Rotation{T}, B::Rotation{T}) where T =
    CompositeRotation([B, A])

Base.:(*)(A::CompositeRotation, B) = rmul!(copy(A), B)
Base.:(*)(A::Rotation, B::CompositeRotation) = lmul!(A, copy(B))

function LinearAlgebra.lmul!(R::CompositeRotation, M)
    for r in R.rotations
        lmul!(r, M)
    end
    M
end

Base.Matrix(R::CompositeRotation{T}, (m,n)) where T =
    lmul!(R, Matrix(one(T)*I, (m,n)))
