# Rainflow
[![Build Status](https://travis-ci.org/dhoegh/Rainflow.jl.svg?branch=master)](https://travis-ci.org/dhoegh/Rainflow.jl) [![Coverage Status](https://img.shields.io/coveralls/dhoegh/Rainflow.jl.svg)](https://coveralls.io/r/dhoegh/Rainflow.jl)
```Julia
using PyPlot
using Rainflow

signal = randn(100) # Gennerates some random data
ext, t = sort_peaks(signal) # Sorts the signal to extremums
plot([0:length(signal)-1], signal)
plot(t,ext,"ro")
cycles = count_cycles(ext,t) # find all the cycles in the data
plot(cycles) # plot of a whole cycle is not plottet correctly, it plots a cylce from the starting point to where the value that defines the range occur.
figure()
plot(cycles[1])
bins = sum_cycles(cycles,10,1) # Sums the cycles together dependant on intervals, here there is 10 range interval and one mean interval
figure()
bar([1:length(bins)],squeeze(bins,2),0.75)

range_intervals = [0, 40, 45, 50, 55, 60, 100.] # There can also be specified user defined intervals
mean_intervals = [0.,100.] # The intevals needs to go from 0-100
bins = sum_cycles(cycles, range_intervals, mean_intervals) # Sums the cycles in the given intervals
figure()
bar([1:length(bins)],squeeze(bins,2),0.75)
```

