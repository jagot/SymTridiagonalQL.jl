module SymTridiagonalQL

using LinearAlgebra
using Random
using Formatting

function horizontal_line(io::IO=stdout; char="━", color=:light_black)
    ncol = displaysize(io)[2]
    printstyled(io, repeat(char, ncol), color=color)
    println(io)
end

include("shifts.jl")
include("triql.jl")

end
