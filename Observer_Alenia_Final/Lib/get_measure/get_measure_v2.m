%% OUTPUT FUNCTION
function [buf_dy_out, yhat] = get_measure_v2(x,j,forward,buf_dy,params)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
%%%%%%%%%%%%%%%

global DynOpt

x_propagate = x;

% check if measure or propagate
if j == 0
    measure_flag = 1;
else
    measure_flag = 0;
end

% start from x0 --> propagate up to x at the jth window
% set j to 0 to measure instantaneously. This is done only if measure_flag
% == 0 namely only during optimisation
if measure_flag == 0
    % get adaptive buffer distance. The zero is added to make the vector of the
    % same length of the buffer as the diff command computes ((i+1)th-ith)
    buf_dist = [max(DynOpt.Y_space(1)-DynOpt.BackTimeIndex,0), diff(DynOpt.Y_space)];
    % get first nonzero element
    [~,pos] = find(buf_dist);
    % check if exist a zero
    zero_flag = pos(1)>1;
    
    if (forward==1)
        % the max has been added to at least save the current state if
        % diff(j) == 0 (no propagation at all)
%         for i=1:max(sum(buf_dist(end-j+1:end)),1)-1
        for i=1:sum(buf_dist(1:pos(j)))-zero_flag
            
            % integration
            [x_propagate, params] = model_propagate_local(DynOpt.BackTimeIndex+i,DynOpt.Ts,x_propagate,params);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update buffer
            x_read = zeros(DynOpt.dim_out,1);
            for k=1:length(params.observed_state)
                x_read(k) = x_propagate(params.observed_state(k));
            end
            % derivatives
            dy = zeros(DynOpt.dim_out,1);
            for k=1:DynOpt.dim_out
                [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,x_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
    else
        %%%%% TO BE DONE %%%%
    end
end

%%%% FIRST OBSERVATION
if measure_flag == 1
    %%%% FIRST OBSERVATION
    x_read = zeros(DynOpt.dim_out,1);
    for i=1:length(params.observed_state)
        x_read(i) = x_propagate(params.observed_state(i));
    end 
    % derivatives
    dy = zeros(DynOpt.dim_out,1);
    for k=1:DynOpt.dim_out
        [buf_dy(k,:), dy(k)] = IterativePseudoDerivative(DynOpt.Ts,x_read(k),DynOpt.c1_derivative,DynOpt.d1_derivative,0,buf_dy(k,:));
    end
end

%%%% 1 OUTPUT %%%%
try 
    if(forward == 1)
        yhat = [x_read, dy];
    else
        yhat = [x_read, -dy];
    end
catch
   disp('ARARMAX');
end

%%%% Buffer output %%%%
buf_dy_out = buf_dy;

end