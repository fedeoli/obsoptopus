%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DynOpt = scale_factor(DynOpt)
    %%%%%%%%%%%%%% reset scale_factor depending on the new buffers %%%%%%%%%%%%
    %%%%% SCALE FACTOR %%%%%
    DynOpt.y_end = 3;
    DynOpt.nJ_nl = 1;
    %%%% bias values %%%%%
    temp_scale = [1;1;0;1];
%     temp_scale = [1;1*1e-1;1*5e-2;5e-1];
    DynOpt.scale_factor_init = 1e0.*temp_scale.*ones(DynOpt.y_end+DynOpt.nJ_nl,9);
    DynOpt.scale_factor_init(:,DynOpt.dim_out+1:end) = 0;
    % memory factor
    % struct.lambda = [0.8; 1; 0.7; 1; 1];
    DynOpt.lambda = [1; 1; 1; 1];
    DynOpt.y_weight = [1*ones(1,3) 1*ones(1,3*DynOpt.nMagneto)];
    DynOpt.scale_factor = zeros(DynOpt.w,DynOpt.y_end+DynOpt.nJ_nl,9);
    for z = 1:DynOpt.w
        DynOpt.scale_factor(DynOpt.w+1-z,:,:) = DynOpt.scale_factor_init.*DynOpt.lambda.^(z-1);
    end
    for z = 1:DynOpt.dim_out
        DynOpt.scale_factor(:,:,z) = DynOpt.scale_factor(:,:,z)*DynOpt.y_weight(z);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%