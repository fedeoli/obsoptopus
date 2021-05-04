%% OUTPUT FUNCTION
function yhat = get_measure_RL_v1(x,params)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
%%%%%%%%%%%%%%%

global DynOpt

x_propagate = x;

%%%% FIRST OBSERVATION
x_read = zeros(DynOpt.dim_out,1);
for i=1:length(params.observed_state)
    x_read(i) = x_propagate(params.observed_state(i));
end 

% get magnetometers
if DynOpt.nMagneto > 0
    [MagTemp_body1, MagTemp_body2] = Magnetometers_GetMeasures_1(DynOpt,params.SatellitesCoordinates(1:3)',x_propagate(1:4)',DynOpt.RPYbetweenMagSensors');
end

if DynOpt.nMagneto == 1
    Mag = MagTemp_body1;
elseif DynOpt.nMagneto == 2
    Mag = [MagTemp_body1; MagTemp_body2];
else
    Mag = [];
end

% output computation
x_read(4:end) = Mag;
y_read = x_read;

%%%% 1 OUTPUT %%%%
try 
   yhat = y_read;
catch
   disp('ARARMAX');
end

end