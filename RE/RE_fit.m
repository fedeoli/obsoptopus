%% fit data
function [signal_interp, time_interp] = RE_fit(time,signal,plot_flag,undersample)

    % interpolate data
    Ts = undersample*max(diff(time));
    time_interp = time(1):Ts:time(end);
    signal_interp = interp1(time,signal,time_interp,'spline');
    
    if plot_flag
       figure(1)
       hold on
       grid on
       plot(time,signal,'o');
       plot(time_interp,signal_interp,'.');
       legend('s','s_int')
    end
end