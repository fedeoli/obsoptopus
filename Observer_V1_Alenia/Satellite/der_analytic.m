%% script to compute the analytical deriative in the main program
%%%% ESTIMATED measurements
% measures
function dy = der_analytic(time_instant,state)

    global DynOpt params
    
    dy = zeros(DynOpt.dim_out,1);
    temp_pos = params.SatellitesCoordinates;
    set_input(time_instant,state,temp_pos);
    dw_temp = AttitudeDynamics_V2_2(state, params);
    dy(1:3) = dw_temp(5:7)*DynOpt.Ts;
end