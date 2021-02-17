%%%%%%%%%%%%%%%%% FIRST MEASURE UPDATE %%%%%%%%
% manage measurements
% set the derivative buffer as before the optimisation process (multiple f computation)
back_time = DynOpt.BackTimeIndex;
if (back_time) >= DynOpt.d1_derivative
    DynOpt.buf_dyhat_temp = DynOpt.Yhat_full_story(:,back_time-(DynOpt.d1_derivative-1):back_time);
else
    init_pos = DynOpt.d1_derivative-back_time;
    DynOpt.buf_dyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,back_time-(back_time-1):back_time)];
end
if (back_time) >= DynOpt.d1_dderivative
    DynOpt.buf_ddyhat_temp = DynOpt.dYhat_full_story(:,back_time-(DynOpt.d1_dderivative-1):back_time);
else
    init_pos = DynOpt.d1_dderivative-back_time;
    DynOpt.buf_ddyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.dYhat_full_story(:,back_time-(back_time-1):back_time)];
end
%%%% ESTIMATED measurements
% measures       
% NB: the output storage has to be done in
% back_time+1 as the propagation has been
% performed 
[DynOpt.buf_dyhat_temp, DynOpt.buf_ddyhat_temp, Yhat] = DynOpt.get_measure(x_propagate,0,measure_forward,DynOpt.buf_dyhat_temp,DynOpt.buf_ddyhat_temp);
DynOpt.Yhat_full_story(:,back_time+1) = Yhat(:,1);
DynOpt.dYhat_full_story(:,back_time+1) = Yhat(:,2);
DynOpt.ddYhat_full_story(:,back_time+1) = Yhat(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% PROPAGATION %%%%%%%%%%%%%%%%%%%%%%%
for j =1:DynOpt.WindowSamples-1     
    % back time
    back_time = DynOpt.BackTimeIndex+j;

    % set input
    set_input(back_time);

    % integrate
    x_propagate = PlantJumpMap_general_notime_params(x_propagate, DynOpt.model, 1, params);
    DynOpt.OptXstory(:,back_time) = x_propagate;

    % manage measurements
    % set the derivative buffer as before the optimisation process (multiple f computation)
    if (back_time) >= DynOpt.d1_derivative
        DynOpt.buf_dyhat_temp = DynOpt.Yhat_full_story(:,back_time-(DynOpt.d1_derivative-1):back_time);
    else
        init_pos = DynOpt.d1_derivative-back_time;
        DynOpt.buf_dyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.Yhat_full_story(:,back_time-(back_time-1):back_time)];
    end
    if (back_time) >= DynOpt.d1_dderivative
        DynOpt.buf_ddyhat_temp = DynOpt.dYhat_full_story(:,back_time-(DynOpt.d1_dderivative-1):back_time);
    else
        init_pos = DynOpt.d1_dderivative-back_time;
        DynOpt.buf_ddyhat_temp = [zeros(DynOpt.dim_out,init_pos), DynOpt.dYhat_full_story(:,back_time-(back_time-1):back_time)];
    end

    %%%% ESTIMATED measurements
    % measures       
    % NB: the output storage has to be done in
    % back_time+1 as the propagation has been
    % performed 
    [DynOpt.buf_dyhat_temp, DynOpt.buf_ddyhat_temp, Yhat] = DynOpt.get_measure(x_propagate,0,measure_forward,DynOpt.buf_dyhat_temp,DynOpt.buf_ddyhat_temp);
    DynOpt.Yhat_full_story(:,back_time+1) = Yhat(:,1);
    DynOpt.dYhat_full_story(:,back_time+1) = Yhat(:,2);
    DynOpt.ddYhat_full_story(:,back_time+1) = Yhat(:,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%