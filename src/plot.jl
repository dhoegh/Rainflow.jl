import PyPlot.plot

function plot(cycle::Cycle)
    time = linspace(cycle.t_s, cycle.t_e, 25)
    amplitude = abs(cycle.v_s-cycle.v_e)/2
    dt = cycle.t_s-cycle.t_e
    if cycle.count==1
        plot(time, cycle.mean.+sign(cycle.v_s-cycle.v_e).*amplitude.*cos.(2π/dt.*(time-cycle.t_s)))
    else
        plot(time, cycle.mean.+sign(cycle.v_s-cycle.v_e).*amplitude.*cos.(π/dt.*(time-cycle.t_s)))
    end
end
