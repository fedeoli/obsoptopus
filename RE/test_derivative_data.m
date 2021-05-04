%% test 
function d1_signal = test_derivative_data(time,signal,plot_flag,d,c)

    % simulation time
    Ts = max(diff(time));
    Niter = length(time);

    % filter setup
    buf_d1 = zeros(1,d);
    
    % output signal
    d1_signal = zeros(1,Niter);

    % compute derivative
    for i=1:Niter
       [buf_d1, dy] = IterativePseudoDerivative(Ts,signal(i),c,d,0,buf_d1);
       d1_signal(i) = dy;
    end

    if plot_flag
        figure(1)
        hold on
        grid on
        ax(1) = subplot(2,1,1);
        plot(time,signal)
        legend('s')
        ax(1) = subplot(2,1,2);
        plot(time,d1_signal)
        legend('ds')
        linkaxes(ax,'x');
    end

end