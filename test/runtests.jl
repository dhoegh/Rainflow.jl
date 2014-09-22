using Rainflow
using Base.Test

srand(1234)
signal = randn(100)
ext, t = sort_peaks(signal)
result = find_cycles(ext,t)
bins = calc_sum(result,13,1)
@test squeeze(bins,2) == [4.0, 6.0, 7.5, 3.0, 2.5, 5.0, 0.0, 2.5, 1.0, 0.0, 0.5, 0.0, 1.0]
