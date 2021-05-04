function X = rk4_V1_1_decentralized(f, tspan, x0, params,ith_satellite)

N = length(x0);             % number of state's components
M = length(tspan);         % numer of time steps
dt = tspan(2) - tspan(1);   % time step

% Matrices allocation
X = zeros(N,M);
X(:,1) = x0;
K1 = zeros(N,M);

for i = 1:M-1
    
    % State and time at ti
    x = X(:,i);
    
    % Runge Kutta 4
    K1(:,i) = f(x, params,ith_satellite,i);
    K2 = f(x + K1(:,i)*dt/2, params,ith_satellite,i);
    K3 = f(x + K2*dt/2, params,ith_satellite,i);
    K4 = f(x + K3*dt, params,ith_satellite,i);
    
    % Solution at ti+1
    X(:,i+1) = x + (dt/6)*(K1(:,i) + 2*K2 + 2*K3 + K4);
    
end


end