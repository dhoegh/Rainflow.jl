module Rainflow

export sort_peaks, find_boundary_vals, count_cycles, sum_cycles, rainflow, rainflow!

# Interval{T} = Union{Vector{T}, StepRangeLen{T}}
struct Cycle  # This is the information stored for each cycle found
    count  ::Float16  # only 1 or 0.5 
    range  ::Float32
    mean   ::Float32
    Rvalue ::Float32   # ratio min(v_s, v_e) / max(v_s, v_e)
    v_s    ::Float32   # value start
    t_s    ::Float32   # time start
    v_e    ::Float32   # value end
    t_e    ::Float32   # time end
    dt     ::Float32   # t_e - t_s
end
end

mutable struct Cycles_bounds #
    min_mean  ::Float64
    max_mean  ::Float64
    max_range ::Float64
    min_R     ::Float64
    max_R     ::Float64
end

Base.show(io::IO,x::Cycle) = print(io, "Cycle: count=", x.count, ", range=",x.range, ", mean=",x.mean, ", R=", x.Rvalue, ", dt=", x.dt)

Base.show(io::IO,x::Cycles_bounds) = print(io, "Cycles_bounds : min mean value=", x.min_mean, ", max mean value=", x.max_mean, ", max range=",x.max_range, ", min R=", x.min_R, ", max R=",x.max_R)

""" This function sort out points where the slope is changing sign"""
function sort_peaks(signal::AbstractVector, dt=collect(1.:length(signal)))
    slope = diff(signal)
    # Determines if the point is local extremum
    is_extremum = vcat(true, (slope[1:end-1].*slope[2:end]).<=0., true)
    return signal[is_extremum] , dt[is_extremum]
end

function cycle(count, v_s, t_s, v_e, t_e)
    Cycle(count, abs(v_s-v_e), (v_s+v_e)/2, min(v_s,v_e)/max(v_s,v_e), v_s, t_s, v_e, t_e, t_e-t_s)
end
end

""" Count the cycles from the data """
function count_cycles(peaks::Vector{T}, t::Vector{S}) where {T, S <: Real}
    list = copy(peaks) # Makes a copy because there is going to be sorted in the vectors
    time = copy(t)

    currentindex = 1
    nextindex    = 2
    cycles       = Vector{Cycle}()
    sizehint!(cycles, length(list)) # This reduces the memory consumption a bit
    
    @inbounds begin
        while length(list) > (currentindex+1)
                currentvalue = abs(list[currentindex+1] - list[currentindex])
                nextvalue    = abs(list[nextindex+1]    - list[nextindex])
            if nextvalue > currentvalue
                if currentindex == 1 # This case counts a half and cycle deletes the poit that is counted
                    push!(cycles, cycle(0.5, list[currentindex], time[currentindex], list[currentindex+1], time[currentindex+1]))
                    popfirst!(list) # Removes the first entrance in ext and time
                    popfirst!(time)
                else # This case counts one cycle and deletes the point that is counted
                    push!(cycles, cycle(1. , list[currentindex], time[currentindex], list[currentindex+1], time[currentindex+1]))
                    deleteat!(list, currentindex:(currentindex+1)) # Removes the i and i+1 entrance in ext and time
                    deleteat!(time, currentindex:(currentindex+1))
                end
                currentindex  = 1
                nextindex     = 2
            else
                currentindex += 1
                nextindex    += 1
            end
        end
        for currentindex = 1:length(list)-1 # This counts the rest of the points that have not been counted as a half cycle
            push!(cycles, cycle(0.5, list[currentindex], time[currentindex], list[currentindex+1], time[currentindex+1]))
        end
    end
    return cycles
end

""" Find the minimum and maximum mean value and maximum range from a vector of cycles"""
function find_boundary_vals(cycles::Vector{Cycle})
    bounds = Cycles_bounds(Inf, -Inf, -Inf, Inf, -Inf)
    for cycle in cycles
        cycle.mean   > bounds.max_mean  && setfield!(bounds, :max_mean,  cycle.mean)
        cycle.mean   < bounds.min_mean  && setfield!(bounds, :min_mean,  cycle.mean)
        cycle.range  > bounds.max_range && setfield!(bounds, :max_range, cycle.range)
        cycle.Rvalue < bounds.min_R     && setfield!(bounds, :min_R,     cycle.Rvalue)
        cycle.Rvalue > bounds.max_R     && setfield!(bounds, :max_R,     cycle.Rvalue)
    end
    return bounds
end

""" Returns the range index the value is belonging in """
function find_range(interval::Vector{T}, value) where {T <: Real}
    for i = 1:length(interval)-1
        if interval[i] <= value < interval[i+1]
            return i
        end
    end
    error("The value where not in range")
end

""" Sums the cycle count given intervals of range_intervals and mean_intervals. The range_intervals and mean_intervals is given in fraction of range size"""
function sum_cycles(cycles::Vector{Cycle}, range_intervals::Interval{T}, mean_intervals::Interval{T}) where {T <: Real}
    bounds = find_boundary_vals(cycles)
    bins   = zeros(length(range_intervals)-1, length(mean_intervals)-1)

    range_in = (range_intervals *  bounds.max_range)/100
    mean_in  = (mean_intervals  * (bounds.max_mean - bounds.min_mean))/100
    mean_in  = mean_in .+ bounds.min_mean

    issorted(mean_intervals)  || error("The array needs to be sorted in raising order")
    issorted(range_intervals) || error("The array needs to be sorted in raising order")

    nr_digits     = 14  # The rounding is performed due to numerical noise in the floats when comparing
    mean_i        = collect(mean_in)
    range_i       = collect(range_in)
    #ensure the cycles is in the intervals by adding a small error to the end values of the interal.
    error_m       = (bounds.max_mean-bounds.min_mean)*1e-14
    mean_i[end]  += error_m
    mean_i[1]    -= error_m
    error_r       =  bounds.max_range*1e-14
    range_i[end] += error_r
    range_i[1]   -= error_r
    #show(mean_intervals)
    for cycle in cycles
        i = find_range(range_i, cycle.range)
        j = find_range(mean_i, cycle.mean)
        bins[i,j] += cycle.count
    end

rainflow(x::AbstractVector, args...) = rainflow(1:1:length(x), x, args...)
function rainflow(t::AbstractVector, x::AbstractVector, args...)
    y, t1  = sort_peaks(x, t)
    res    = count_cycles(y, t1)
    bounds = find_boundary_vals(res)
    range_i, mean_i = set_rangemx(bounds, args...)
    bins   = zeros(length(range_i)-1, length(mean_i)-1)
    sum_cycles!(bins, res, range_i, mean_i)
    return bins
end

function rainflow!(bins::AbstractMatrix, t::AbstractVector, x::AbstractVector, range_i::AbstractVector, mean_i::AbstractVector, args...)
    y, t1  = sort_peaks(x, t)
    res    = count_cycles(y, t1)
    sum_cycles!(bins, res, range_i, mean_i, args...)
end

function rainflow(t, x, args...)
    y, t1 = sort_peaks(x, t)
    res   = count_cycles(y, t1)
    bins  = sum_cycles(res, args...)
end

end # module