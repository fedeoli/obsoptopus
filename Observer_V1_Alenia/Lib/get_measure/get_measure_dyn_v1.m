%% OUTPUT FUNCTION
function yhat = get_measure_dyn_v1(x,j,forward,params)

% Output function to retrieve the estimated output from the estimated state
% x and the input u (linear case yhat = C*x + D*u)
%%%%%%%%%%%%%%%

global DynOpt

x_propagate = x;

% check if measure or propagate
measure_flag = 0;

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
    % set n_iter
    n_iter = sum(buf_dist(1:pos(j)))-zero_flag;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (forward==1)
        % the max has been added to at least save the current state if
        % diff(j) == 0 (no propagation at all)
%         for i=1:max(sum(buf_dist(end-j+1:end)),1)-1
        for i=1:n_iter
            
            % integration
            [x_propagate, params] = model_propagate_local(DynOpt.BackTimeIndex+i,DynOpt.Ts,x_propagate,params);
            
            % dynamics
            xdot = DynOpt.model(x_propagate,params);

        end
    else
        %%%%% TO BE DONE %%%%
    end
end

%%%% 1 OUTPUT %%%%
try 
    if(forward == 1)
        yhat = xdot;
    else
        yhat = xdot;
    end
catch
   disp('ARARMAX');
end

end