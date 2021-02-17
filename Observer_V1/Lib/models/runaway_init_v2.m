%% PLANT model and data
% model simulation
% plant data
function runaway_init_v2(struct)
    % global vars
    global DynOpt params 
    
    % get measurement data
    if (DynOpt.simulationModel == 0) || (DynOpt.check == 1)
        disp('Get data')
        DynOpt.Tstart = struct.T0;
        DynOpt.Tend = struct.Tend;
        DynOpt.sample_time = struct.sample_time;
        DynOpt.data = struct.dati;

        [~,cols] = size(DynOpt.data);

        [~,pos_0] = min(abs(DynOpt.data(:,1)-struct.T0));
        [~,pos_end] = min(abs(DynOpt.data(:,1)-struct.Tend));

        time_sample = pos_0:DynOpt.sample_time:pos_end;
        values = struct.dati(time_sample,2);

        DynOpt.dataStory = zeros(cols-1,length(values));
        DynOpt.dataStory(1,:) = values;
        
        % oversample
        DynOpt.Ts = 1e-1*DynOpt.Ts;
        XQ = DynOpt.time(1):DynOpt.Ts:DynOpt.time(end);
        YQ = pchip(DynOpt.time,DynOpt.dataStory,XQ);
        
        DynOpt.dataStory_old = DynOpt.dataStory;
        DynOpt.time_old = DynOpt.time;
        
        DynOpt.dataStory = YQ;
        DynOpt.time = XQ;
        DynOpt.Niter = length(DynOpt.time);
        
        if DynOpt.check == 0
            DynOpt.time = DynOpt.data(time_sample);
        end
    end

    % electric charge
    params.Q = 190;

    % spontaneous emission
    params.S = 0;
    
    % eps_coef
    params.eps_coef = 5.65e2;
    
    % ringing
    params.wq = 30;
    params.chi = 20;
    params.c1 = 1.1e3;
    params.offset = 50;
    params.chi_offset = 20;

    % parameters
    if struct.identify == 1
        params.gamma = 0.75;
        params.gamma1 = 60;
        params.ni = 6e1;
        params.Wt = 1e1;
    else
        params.gamma = struct.gamma;
        params.gamma1 = struct.gamma1;
        params.ni = 0.5;
        params.Wt = struct.Wt;
    end


    % initial condition
    if DynOpt.simulationModel == 1
        params.T0 = 150;
        if DynOpt.check == 1
           params.W0 = DynOpt.dataStory(1,1);
%            params.W0 = 0.001;
        else
           params.W0 = 0.001; 
        end
        params.x10 = 0;
        params.x20 = 0;
    else
        params.T0 = 0;
        params.W0 = DynOpt.dataStory(1,1);
        params.x10 = 0;
        params.x20 = 0;
    end
    
    %%%%%% STEADY STATE %%%%%%
%     delta = (2*params.ni*params.Wt-params.gamma*params.Q+2*params.gamma1*params.Wt)^2+8*params.ni*params.gamma*params.Q*params.Wt;
%     params.Wss1 = (-(2*params.ni*params.Wt-params.gamma*params.Q+2*params.gamma1*params.Wt)+sqrt(delta))/(4*params.ni);
%     params.Tss1 = (params.Q-params.S)/(2*params.Wss1);
%     params.Wss2 = (-(2*params.ni*params.Wt-params.gamma*params.Q+2*params.gamma1*params.Wt)-sqrt(delta))/(4*params.ni);
%     params.Tss2 = (params.Q-params.S)/(2*params.Wss2);

    %%%%%%%%% WORKSPACE %%%%%%%%%%
    % params.T0 = 9.255;
    % params.W0 = 0.02176;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     DynOpt.param_estimate = [params.gamma, params.gamma1, params.Q, params.eps_coef, params.wq, params.chi];%, params.ni, params.Wt, params.Q, params.S];
    DynOpt.param_estimate = [params.Q,params.eps_coef,params.wq, params.chi];%, params.ni, params.Wt, params.Q, params.S];

    params.StateDim = 4;

    if DynOpt.simulationModel == 1
        params.observed_state = [2,3];
    else
        params.observed_state = [2,3];
    end

    params.input_flag = 0;

    if params.input_flag == 0
        params.U = zeros(1,DynOpt.Niter);
    else
        params.U = zeros(1,DynOpt.Niter);
    end

    DynOpt.params = params;
end