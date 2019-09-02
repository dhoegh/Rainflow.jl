# Rainflow-counting algorithm
[![Build Status](https://travis-ci.org/dhoegh/Rainflow.jl.svg?branch=master)](https://travis-ci.org/dhoegh/Rainflow.jl) [![Coverage Status](https://img.shields.io/coveralls/dhoegh/Rainflow.jl.svg)](https://coveralls.io/r/dhoegh/Rainflow.jl)

This package contains an implementation of the ASTM E1049-85 rainflow-counting algorithem and it can be used to asses fatigue properties of a structure by reducing a stress spectrum down to simple stress reversals that can be handeled by Palmgren-Miner linear damage rule.
# Installation

The package can be installed by:
```julia
Pkg.add("Rainflow")
```
And a plotting package to visualize the results:
```julia
Pkg.add("Plots")
```
# Use
The use of the package can be seen in the following example:
```Julia
using Rainflow
using Plots

signal = 10*randn(100); # Gennerates some random data
extremum, t = sort_peaks(signal) # Sorts the signal to local extrema's, could optional take a time vector
plot(signal)
scatter!(t,extremum) # plots extrema's
cycles = count_cycles(extremum, t) # find all the cycles in the data

bins = sum_cycles(cycles, 10, 1) # Sums the cycles together dependant on intervals, here there is 10 range interval and one mean interval
bar(bins)

range_intervals = [0, 40, 45, 50, 55, 60, 100.] # There can also be specified user defined intervals
mean_intervals = [0.,100.] # The intevals needs to go from 0-100
bins = sum_cycles(cycles, range_intervals, mean_intervals) # Sums the cycles in the given intervals

# Plotting the results with the user defined intervals
boundaries = find_boundary_vals(cycles); # Find the min max values to scale or interpolate the calculated values
range_values = range_intervals/100*boundaries.max_range # Scale the range values
# The bar plot in julia only takes integer array as tick marks
range_values = collect(Iterators.drop(range_values,1)) # Drop the interval with 0 amplitude
range_values = map(x->round.(x), range_values) # round off the values so they can be converted to Int
range_values = convert(Array{Int64,1},range_values); # Convert the array to Int
bar(range_values, bins, xticks = range_values)
```

## Performance note
The algorithem and functions has been tested on data set of 1e8 numbers, and the sorting and sum the cycles into bins can be performed in less than 20s on a laptop.

Don't hesitate to file an issue or pull-request to improve the package.
