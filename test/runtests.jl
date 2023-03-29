import Pkg
Pkg.activate(normpath(@__DIR__, ".."))

using Rainflow
using Test

t = [0, 1, 3, 4, 5, 6, 8, 10, 13, 15]
x = [-2, 1, -3, 5, -1, 3, -4, 4, -2, 6]
ext, t = sort_peaks(x, t)
result = count_cycles(ext, t)
bins   = sum_cycles(result, 10, 1)
bins1  = sum_cycles(result, collect(range(0, 100, length=10+1)), collect(range(0, 100, length=1+1)))
bins2  = sum_cycles(result, range(0, 100, length=10+1), range(0, 100, length=1+1))
bins3  = rainflow(t, x)

@testset "Rainflow" begin
    @test bins == bins1 == bins2 == bins3
    @test getfield.(result, :count) == [0.5, 0.5, 1, 0.5, 1, 0.5, 0.5]
    @test dropdims(bins, dims = 2)  == [0.0, 0.0, 0.0, 0.5, 1.5, 0.0, 1.0, 0.0, 0.5, 1.0]
end