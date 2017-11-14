module Rainflow
import Base.show

export sort_peaks, find_boundary_vals, count_cycles, sum_cycles

""" This function sort out points where the slope is changing sign."""
function sort_peaks(signal::AbstractArray{Float64,1}, dt=collect(0.:length(signal)-1.))
    slope = diff(signal)
    # Determines if the point is local extremum
    is_extremum = vcat(true, (slope[1:end-1].*slope[2:end]).<=0., true)
    return signal[is_extremum] , dt[is_extremum]
end

struct Cycle  # This is the information stored for each cycle found
    count::Float64
    range::Float64
    mean::Float64  # value
    v_s::Float64   # value start
    t_s::Float64   # time start
    v_e::Float64   # value end
    t_e::Float64   # time end
end

show(io::IO,x::Cycle) = print(io, "Cycle: count=", x.count, ", range=",x.range, ", mean=",x.mean)

function cycle(count::Float64, v_s::Float64, t_s::Float64, v_e::Float64, t_e::Float64)
    Cycle(count, abs(v_s-v_e), (v_s+v_e)/2, v_s, t_s, v_e, t_e)
end

""" Count the cycles from the data. """
function count_cycles(ext_in::Array{Float64,1},t::Array{Float64,1})
    ext = copy(ext_in) # Makes a copy because there is going to be sorted in the vectors
    time = copy(t)
    i = 1
    j = 2
    cycles = Cycle[]
    sizehint!(cycles, length(ext)) # This reduces the memory consumption a bit
    @inbounds begin
    while length(ext)>(i+1)
        Y = abs(ext[i+1]-ext[i])
        X = abs(ext[j+1]-ext[j])
        if X>=Y
            if i == 1 # This case counts a half and cycle deletes the poit that is counted
                #println("Half cycle $(ext[i]), $(ext[i+1])")
                push!(cycles,cycle(0.5 ,ext[i], time[i], ext[i+1],time[i+1]))
                shift!(ext) # Removes the first entrance in ext and time
                shift!(time)
            else # This case counts one cycle and deletes the poit that is counted
                #println("One cycle $(ext[i]), $(ext[i+1])")
                push!(cycles,cycle(1. ,ext[i], time[i], ext[i+1], time[i+1]))  
                # Removes the i and i+1 entrance in ext and time
                my_deleteat!(ext, i:(i+1))
                my_deleteat!(time, i:(i+1))
            end
            i = 1
            j = 2
        else
            i += 1
            j += 1
        end
    end
    for i=1:length(ext)-1 # This counts the rest of the points that have not been counted as a half cycle
        #println("Half cycle $(ext[i]), $(ext[i+1])")
        push!(cycles,cycle(0.5 ,ext[i], time[i], ext[i+1],time[i+1]))
    end
    end
    return cycles
end

mutable struct Cycles_bounds #
    min_mean::Float64
    max_mean::Float64
    max_range::Float64
end

show(io::IO,x::Cycles_bounds) = print(io, "Cycles_bounds : min mean value=", x.min_mean, ", max mean value=", x.max_mean, ", max range=",x.max_range)

""" Find the minimum and maximum mean value and maximum range from a vector of cycles. """
function find_boundary_vals(cycles::Array{Cycle,1})
    bounds = Cycles_bounds(Inf, -Inf, -Inf)
    for cycle in cycles
        cycle.mean > bounds.max_mean && setfield!(bounds, :max_mean, cycle.mean)
        cycle.mean < bounds.min_mean && setfield!(bounds, :min_mean, cycle.mean)
        cycle.range > bounds.max_range && setfield!(bounds, :max_range, cycle.range)
    end
    return bounds
end

""" Returns the range index the value is belonging in """
function find_range{T<:Real}(interval::Array{T,1},value)
    for i=1:length(interval)-1
        if interval[i] <= value < interval[i+1]
            return i
        end
    end
    error("The value where not in range")
end

    
Interval{T} = Union{Array{T,1}, StepRangeLen{T}}

""" Sums the cycle count given intervals of range_intervals and mean_intervals. The range_intervals and mean_intervals is given in fraction of range size"""
function sum_cycles{T<:Real}(cycles::Array{Cycle,1}, range_intervals::Interval{T}, mean_intervals::Interval{T})
    bounds = find_boundary_vals(cycles)
    bins = zeros(length(range_intervals)-1, length(mean_intervals)-1)
    range_in = (range_intervals*bounds.max_range)/100
    mean_in = (mean_intervals*(bounds.max_mean-bounds.min_mean))/100
    mean_in += bounds.min_mean
    issorted(mean_intervals) || error("The array needs to be sorted in raising order")
    issorted(range_intervals) || error("The array needs to be sorted in raising order")
    nr_digits = 14  # The rounding is performed due to numerical noise in the floats when comparing
    mean_i = collect(mean_in)
    range_i = collect(range_in)
    #ensure the cycles is in the intervals by adding a small error to the end values of the interal.
    error_m = (bounds.max_mean-bounds.min_mean)*1e-14
    mean_i[end]+=error_m
    mean_i[1]-=error_m
    error_r = bounds.max_range*1e-14
    range_i[end]+=error_r
    range_i[1]-=error_r
    #show(mean_intervals)
    for cycle in cycles
        i = find_range(range_i, cycle.range)
        j = find_range(mean_i, cycle.mean)
        bins[i,j] += cycle.count
    end
    return bins
end

function sum_cycles(cycles::Array{Cycle,1}, nr_ranges::Int=10, nr_means::Int=1)
    range_intervals = linspace(0,100,nr_ranges+1)
    mean_intervals = linspace(0,100,nr_means+1)
    sum_cycles(cycles, range_intervals, mean_intervals)
end

# This is necersary because of https://github.com/JuliaLang/julia/issues/24494
if VERSION<v"0.6.2" && Base.is_windows()
    function my_deleteat!(a::Vector, r::UnitRange{<:Integer})
        n = length(a)
        isempty(r) || _deleteat_beg!(a, first(r), length(r))
        return a
    end
    function _deleteat_beg!(a::Vector, i::Integer, delta::Integer)
        if i > 1
            ccall(:memmove, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Csize_t),
                  pointer(a, 1+delta), pointer(a, 1), (i-1)*Base.elsize(a))
        end
        ccall(:jl_array_del_beg, Void, (Any, UInt), a, delta)
        return a
    end
else
    const my_deleteat! = deleteat!
end

try
    include("plot.jl")
catch
    println("""For added plotting features do: Pkg.add("PyPlot")""")
end

end # module
