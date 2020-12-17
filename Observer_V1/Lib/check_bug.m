%% script to check bug (see sensitivity issue)
% check bug 
global DynOpt
if strcmp(DynOpt.opt_method,'sensitivity')
        windowsamples = floor(DynOpt.WindowSamples/(DynOpt.Nts+1));   
        checkpoint = max(0,floor(DynOpt.ActualTimeIndex/(DynOpt.Nts+1))-windowsamples);
        BackTimeIndex_back = (DynOpt.ActualTimeIndex)-(DynOpt.w)*DynOpt.Nts; 
        BackTimeIndex = (DynOpt.ActualTimeIndex)-(DynOpt.w-1)*DynOpt.Nts; 
        if checkpoint == DynOpt.checkpoint_back & BackTimeIndex_back > 0 & DynOpt.ActualTimeIndex >= max(1,DynOpt.WindowSamples) & mod(DynOpt.ActualTimeIndex,DynOpt.Nts) == 1
           DynOpt.back_flag(end+1) = 1;
           
           % counter reset
           DynOpt.back_flag_counter = 0;
           
           % hold reset
           DynOpt.back_flag_hold = 0;
           
           % shift
           DynOpt.shift_index(end+1) = BackTimeIndex;
           DynOpt.n_shift = length(DynOpt.shift_index)-1;
           
           % state copy
           DynOpt.OptXstory(:,BackTimeIndex:BackTimeIndex+DynOpt.WindowSamples) = DynOpt.OptXstory(:,BackTimeIndex_back:BackTimeIndex_back+DynOpt.WindowSamples);
           DynOpt.Xstory(:,BackTimeIndex:BackTimeIndex+DynOpt.WindowSamples) = DynOpt.Xstory(:,BackTimeIndex_back:BackTimeIndex_back+DynOpt.WindowSamples);
           
           % measure copy
           Y_noise = DynOpt.Y(:,end);
           DynOpt.Y(:,1:end-1) = DynOpt.Y(:,2:end);
           
           x_temp = DynOpt.OptXstory(:,BackTimeIndex+(DynOpt.w-1)*DynOpt.Nts);
           for j=1:DynOpt.Nts
              x_temp = PlantJumpMap_general_notime_params(x_temp,DynOpt.model,DynOpt.ForwardOptimization,params);
           end
           DynOpt.Y(:,end) = x_temp(DynOpt.params.observed_state);
           
           % derivative buffer restore    
           startpoint = 1;
           DynOpt.buf_dy = [zeros(1,startpoint), DynOpt.buf_dy(startpoint+1:end)];
           DynOpt.buf_dyhat = [zeros(1,startpoint), DynOpt.buf_dyhat(startpoint+1:end)];
        else
           DynOpt.back_flag(end+1) = 0;
        end
        
        if mod(DynOpt.ActualTimeIndex,DynOpt.Nts) == 1
            DynOpt.checkpoint_back = checkpoint;
        end
        
        if mod(DynOpt.ActualTimeIndex,DynOpt.Nts) == 1 & length(DynOpt.back_flag) > DynOpt.Nts
            DynOpt.back_flag_old = DynOpt.back_flag(end-DynOpt.Nts);
            
            if DynOpt.back_flag(end) == 1
                DynOpt.back_flag_hold = 1;         
            end
        end
        
        if DynOpt.back_flag_hold == 1 & mod(DynOpt.ActualTimeIndex,DynOpt.Nts) == 1
           DynOpt.back_flag_counter = DynOpt.back_flag_counter+1; 
        end
else
        DynOpt.back_flag = 0;
        DynOpt.n_shift = 0;
        DynOpt.back_flag_counter = 0;
        DynOpt.back_flag_hold = 0;
end