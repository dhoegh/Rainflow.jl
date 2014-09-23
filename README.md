# Rainflow
```Julia
using PyPlot
using Rainflow

srand(1234)
signal = randn(100)
ext, t = sort_peaks(signal)
plot([0:length(signal)-1], signal)
plot(t,ext,"ro")
cycles = find_cycles(ext,t)
plot(cycles)
figure()
plot(cycles[1])
stats = check_max(cycles)
bins = calc_sum(cycles,13,1)
figure()
bar([1:13],squeeze(bins,2),0.75)
xlim(1,14)

L = [0, 40, 45, 50, 55, 60, 100.]
b = calc_sum(cycles,L,[0.,100.])
figure()
bar([1:length(b)],squeeze(b,2),0.75)
xlim(1,length(b))
squeeze(b,2) == [24, 2.0, 2.0, 1.5, 1.0, 2.5]
```
[![Build Status](https://travis-ci.org/dhoegh/Rainflow.jl.svg?branch=master)](https://travis-ci.org/dhoegh/Rainflow.jl)

[![Coverage Status](https://img.shields.io/coveralls/dhoegh/Rainflow.jl.svg)]
