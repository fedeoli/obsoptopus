%% interpolate
plot_flag = 0;
undersample = 0.02;
[Data.signal_interp, Data.time_interp] = RE_fit(Data.time,Data.neu213,plot_flag,undersample);

if plot_flag
    pause()
    close all
end

%% filter
if 0
    [Data.butter.A,Data.butter.B] = butter(1,0.999);
    Data.signal_interp_filter = filter(Data.butter.B,Data.butter.A,Data.signal_interp);

    % plot filtering 
    plot(Data.time_interp,Data.signal_interp,'o');
    hold on
    plot(Data.time_interp,Data.signal_interp_filter,'--');
end

%% derivative
d = 10;
c = 4;
Data.signal_interp_der = test_derivative_data(Data.time_interp,Data.signal_interp,plot_flag,d,c);
if plot_flag
    pause()
    close all
end

%% data setup for curve fitting
Data.n = Data.signal_interp;
Data.dn = Data.signal_interp_der;
[Data.vloop_interp, ~] = RE_fit(Data.time,Data.vloop,0,undersample);
[Data.density_interp, ~] = RE_fit(Data.time,Data.density,0,undersample);

interval = (d):length(Data.time_interp);
n = Data.n(interval);
dn = Data.dn(interval);
vloop = Data.vloop_interp(interval);
density = Data.density_interp(interval);
wrap_term = (vloop./(2*pi*0.96)).*(3e19./density);

%% nonlinear least squares
version = 'dn';
if strcmp(version,'dn')
    b = transpose(dn);
    A = [transpose(n), transpose(wrap_term)];
    theta = lsqr(A,b);

    %%% test lqr - dn %%%
    dn_test = zeros(length(interval),1);
    a = theta(1);
    b = theta(2);
    for i=1:length(interval)
        dn_test(i) = RE_fun_fit(dn(i),n(i),vloop(i),density(i),a,b,version);
    end
    n_test = cumtrapz(Data.time_interp(interval),dn_test)+n(1);

    figure
    hold on
    grid on
    plot(Data.time_interp(interval),n,'o');
    plot(Data.time_interp(interval),n_test,'.');
    legend('data','fit')
elseif strcmp(version,'n')
    b = transpose(n);
    A = [transpose(dn), transpose(wrap_term)];
    theta = lsqr(A,b);
    
    n_test = zeros(length(interval),1);
    a = theta(1);
    b = -a*theta(2);
    for i=1:length(interval)
        n_test(i) = RE_fun_fit(dn(i),n(i),vloop(i),density(i),a,b,version);
%         n_test(i) = theta(1)*dn(i) + theta(2)*wrap_term(i);
    end

    figure
    hold on
    grid on
    plot(Data.time_interp(interval),n,'o');
    plot(Data.time_interp(interval),n_test,'.');
    legend('data','fit')
end