module Rainflow

export sort_peaks, check_max, find_cycles, calc_sum

if Pkg.installed("PyPlot")
  include("plot.jl")
end

function sort_peaks(signal::AbstractArray{Float64,1}, dt=[0.:length(signal)-1.])
    """ This sort out points where the slope is changing sign"""
    slope = diff(signal)
    # Determines if the point is local extremum
    is_extremum = [true ;(slope[1:end-1].*slope[2:end]).<=0.;true] 
    return signal[is_extremum] , dt[is_extremum]
end

type Cycle  # This is the information stored for each cycle found
    count::Float64
    range::Float64
    mean::Float64  # value
    v_s::Float64   # value start
    t_s::Float64   # time start
    v_e::Float64   # value end
    t_e::Float64   # time end
end

function cycle(count::Float64, v_s::Float64, t_s::Float64, v_e::Float64, t_e::Float64)
    Cycle(count, abs(v_s-v_e), (v_s+v_e)/2, v_s, t_s, v_e, t_e)
end

function find_cycles(ext_in::Array{Float64,1},t::Array{Float64,1})
    ext = copy(ext_in) # Makes a copy because there is going to be sorted in the vectors
    time = copy(t) 
    i = 1
    j = 2
    cycles = Cycle[]
    sizehint(cycles, length(ext)) # This reduces the memory consumption a bit
    @inbounds begin
    while length(ext)>(i+1)
        Y = abs(ext[i+1]-ext[i])
        X = abs(ext[j+1]-ext[j])
        if X>=Y
            if i == 1 # This case covers a half cycle deletes the poit that is counted
                #println("Half cycle $(ext[i]), $(ext[i+1])")
                push!(cycles,cycle(0.5 ,ext[i], time[i], ext[i+1],time[i+1]))
                shift!(ext) # Removes the first entrance in ext and time
                shift!(time)
            else # This case covers a one cycle deletes the poit that is counted
                #println("One cycle $(ext[i]), $(ext[i+1])")
                push!(cycles,cycle(1. ,ext[i], time[i], ext[i+1],time[i+1]))
                splice!(ext,i+1)  # Removes the i and i+1 entrance in ext and time
                splice!(ext,i)
                splice!(time,i+1)
                splice!(time,i)
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
    return cycles
    end
end

type cycle_stats # 
    min_mean::Float64
    max_mean::Float64
    max_range::Float64
end

function check_max(cycles::Array{Cycle,1})
    stats = cycle_stats(Inf, -Inf, -Inf)
    for cycle in cycles
        cycle.mean > stats.max_mean && setfield!(stats, :max_mean, cycle.mean)
        cycle.mean < stats.min_mean && setfield!(stats, :min_mean, cycle.mean)
        cycle.range > stats.max_range && setfield!(stats, :max_range, cycle.range)
    end
    return stats
end

function calc_sum{T<:Real}(cycles::Array{Cycle,1}, range_spacing::Array{T,1}, mean_spacing::Array{T,1})
    stats = check_max(cycles)
    bins = zeros(length(range_spacing)-1, length(mean_spacing)-1)
    range_spacing *= stats.max_range/100
    #show(range_spacing)
    mean_spacing *= (stats.max_mean-stats.min_mean)/100
    mean_spacing += stats.min_mean
    #show(mean_spacing)
    for cycle in cycles
        k=false
        for i=1:length(range_spacing)-1
            if round(range_spacing[i],15)< round(cycle.range,15) <= round(range_spacing[i+1],15)
                for j=1:length(mean_spacing)-1
                    if round(mean_spacing[j],15)<= round(cycle.mean,15) <= round(mean_spacing[j+1],15)
                        bins[i,j] += cycle.count
                        k=true
                    end
                end
            end
        end
        if !k
            error("A cylce is missed in the counting")
        end
    end
    return bins
end

function calc_sum(cycles::Array{Cycle,1}, nr_ranges::Int=13, nr_means::Int=1)
    range_spacing = linspace(0,100,nr_ranges+1)
    mean_spacing = linspace(0,100,nr_means+1)
    calc_sum(cycles, range_spacing, mean_spacing)
end

end # module
