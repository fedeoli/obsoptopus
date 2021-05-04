%% derivative of the flow - Khalil chapter 3
function csi2_dot = satellite_dflow_input_v2(t,csi1,csi2,params)

    global DynOpt

    % define vars
    q = csi1(1:4);
    w = csi1(5:7);
    I = csi1(8:10);
    
    % parameters
    r = params.SatellitesCoordinates(1:3);
    r_norm = norm(r);
    v = params.SatellitesCoordinates(4:6);
    v_norm = norm(v);
    
    % control parameters
    q_ref = params.DesiredAttitude;
    
    % subs
    temp_dF = subs(DynOpt.dF_sym,DynOpt.q_sym,q);
    temp_dF = subs(temp_dF,DynOpt.w_sym,w);
    temp_dF = subs(temp_dF,DynOpt.I_vec_sym,I);
    temp_dF = subs(temp_dF,DynOpt.r_sym,r);
    temp_dF = subs(temp_dF,DynOpt.r_norm_sym,r_norm);
    temp_dF = subs(temp_dF,DynOpt.v_sym,v);
    temp_dF = subs(temp_dF,DynOpt.v_norm_sym,v_norm);
    temp_dF = subs(temp_dF,DynOpt.mie_sym,params.mi);
    temp_dF = subs(temp_dF,DynOpt.qr,q_ref);
    temp_dF = double(temp_dF);
    
    % fill
    temp_dF = [temp_dF, zeros(7,3)];
    temp_dF = [temp_dF; zeros(6,13)];
    dflow = temp_dF;
    
    % evolution
    csi2_dot = dflow*csi2;
end