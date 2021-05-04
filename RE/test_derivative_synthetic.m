%% test 
T0 = 0;
Tend = 10;
Ts = 1e-2;
time = T0:Ts:Tend;
Niter = length(time);

d = 5;
c = 3;
buf_d1 = zeros(1,d);
buf_d2 = zeros(1,d);

d1_signal = zeros(1,Niter);
d2_signal = zeros(1,Niter);

signal = 0*sin(2*time)+1*cos(time);

for i=1:Niter
   [buf_d1, dy] = IterativePseudoDerivative(Ts,signal(i),c,d,0,buf_d1);
   d1_signal(i) = wrapToPi(dy);
   
   [buf_d2, ddy] = IterativePseudoDerivative(Ts,d1_signal(i),c,d,0,buf_d2);
   d2_signal(i) = wrapToPi(ddy);
end

figure(1)
hold on
grid on
plot(time,signal)
plot(time,d1_signal)
plot(time,d2_signal)
legend('s','ds','dds')

