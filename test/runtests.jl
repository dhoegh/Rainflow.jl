using Rainflow
using Base.Test

srand(1234)
signal = randn(100)
ext, t = sort_peaks(signal)
result = count_cycles(ext,t)
bins = sum_cycles(result,13,1)
bins1 = sum_cycles(result, collect(linspace(0, 100, 13+1)), collect(linspace(0, 100, 1+1)))
bins2 = sum_cycles(result, linspace(0, 100, 13+1), linspace(0, 100, 1+1))
@test bins==bins1==bins2
@test squeeze(bins,2) == [4.0, 6.0, 7.5, 3.0, 2.5, 5.0, 0.0, 2.5, 1.0, 0.0, 0.5, 0.0, 1.0]

