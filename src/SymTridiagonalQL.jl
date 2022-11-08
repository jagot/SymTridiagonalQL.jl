module SymTridiagonalQL

using LinearAlgebra
using Random
using Formatting

dotu(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}) =  dot(a, b)

function dotu(a::AbstractVector{<:Complex}, b::AbstractVector{<:Complex})
    N = length(a)
    length(b) == N || throw(DimensionMismatch("Lengths of arguments do not agree"))
    LinearAlgebra.BLAS.dotu(N, a, stride(a, 1),
                            b, stride(b, 1))
end

function horizontal_line(io::IO=stdout; char="━", color=:light_black)
    ncol = displaysize(io)[2]
    printstyled(io, repeat(char, ncol), color=color)
    println(io)
end

include("shifts.jl")
include("rotations.jl")
include("triql.jl")

end
