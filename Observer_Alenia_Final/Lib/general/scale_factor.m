%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DynOpt = scale_factor(DynOpt)
    %%%%%%%%%%%%%% reset scale_factor depending on the new buffers %%%%%%%%%%%%
    %%%%% SCALE FACTOR %%%%%
    %%%% bias values %%%%%
    DynOpt.scale_factor_init = 1e0.*DynOpt.temp_scale.*ones(DynOpt.y_end+DynOpt.nJ_nl,9);
    DynOpt.scale_factor_init(:,DynOpt.dim_out+1:end) = 0;
    % memory factor
    DynOpt.scale_factor = zeros(DynOpt.w,DynOpt.y_end+DynOpt.nJ_nl,9);
    for z = 1:DynOpt.w
        DynOpt.scale_factor(DynOpt.w+1-z,:,:) = DynOpt.scale_factor_init.*DynOpt.lambda.^(z-1);
    end
    for z = 1:DynOpt.dim_out
        DynOpt.scale_factor(:,:,z) = DynOpt.scale_factor(:,:,z)*DynOpt.y_weight(z);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%