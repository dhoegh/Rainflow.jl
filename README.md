# Rainflow-counting algorithm
[![Build Status](https://travis-ci.org/dhoegh/Rainflow.jl.svg?branch=master)](https://travis-ci.org/dhoegh/Rainflow.jl) [![Coverage Status](https://img.shields.io/coveralls/dhoegh/Rainflow.jl.svg)](https://coveralls.io/r/dhoegh/Rainflow.jl)

This package contains an implementation of the ASTM E1049-85 rainflow-counting algorithem and it can be used to asses fatigue properties of a structure by reducing a stress spectrum down to simple stress reversals that can be handeled by Palmgren-Miner linear damage rule.
# Installation

The package can be installed by:
```julia
Pkg.add("Rainflow")
```
If plotting of the cycles is desired PyPlot also needs to be installed and I recommend the use of Continuum Analytics python distribution Anaconda
```julia
Pkg.add("PyPlot")
```
# Use
The use of the package can be seen in the following example:
```Julia
using PyPlot
using Rainflow

signal = randn(100) # Gennerates some random data
ext, t = sort_peaks(signal) # Sorts the signal to local extrema's, could optional take a time vector
plot(collect(0:length(signal)-1), signal)
plot(t,ext,"ro") # plots extrema's
cycles = count_cycles(ext, t) # find all the cycles in the data
plot.(cycles) # plot of a whole cycle is not plottet correctly, it plots a cylce from the starting point to where the value that defines the range occur.
figure()
plot(cycles[1])
bins = sum_cycles(cycles, 10, 1) # Sums the cycles together dependant on intervals, here there is 10 range interval and one mean interval
figure()
bar(collect(1:length(bins)), squeeze(bins,2), 0.75)

range_intervals = [0, 40, 45, 50, 55, 60, 100.] # There can also be specified user defined intervals
mean_intervals = [0.,100.] # The intevals needs to go from 0-100
bins = sum_cycles(cycles, range_intervals, mean_intervals) # Sums the cycles in the given intervals
figure()
bar(collect(1:length(bins)), squeeze(bins,2), 0.75)
```

## Performance note
The algorithem and functions has been tested on data set of 1e8 numbers, and the sorting and sum the cycles into bins can be performed in less than 20s on a laptop.

Don't hesitate to file an issue or pull-request to improve the package.
