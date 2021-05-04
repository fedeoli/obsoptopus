%% symbolic analysys
function [DynOpt_out,params_out] = SymAnalysis_RL_v4

    %%% global vars %%%
    global DynOpt params
    
    %%%%%%%%%%%%%%%%%%%% GENERATE MAPS %%%%%%%%%%%%%%%%%%
    %%%%% Sym analysis %%%%%
    fprintf('Setting Observer parameters\n');
    syms Ixx Iyy Izz;       % System Inertia
    syms wx wy wz;          % System angular velocity
    syms q1 q2 q3 q0;       % Quaternion
    syms rx ry rz r_norm;   % sat-Earth vector
    syms vx vy vz v_norm;   % sat-Earth velocity
    syms mie;               % Earth constant
    syms Bx1 By1 Bz1
    syms Bx2 By2 Bz2
    syms dt1 dt2 dt3
    syms taux tauy tauz
    
    symarray_H = [ Ixx Iyy Izz ...
                 wx wy wz ...
                 q1 q2 q3 q0 ...
                 rx ry rz r_norm ...
                 vx vy vz v_norm ...
                 mie ...
                 Bx1 By1 Bz1 ...
                 Bx2 By2 Bz2 ...
                 dt1 dt2 dt3 ...
                 taux tauy tauz];
             
     symarray_G = [ Ixx Iyy Izz ...
                 wx wy wz ...
                 q1 q2 q3 q0 ...
                 rx ry rz r_norm ...
                 vx vy vz v_norm ...
                 mie, ...
                 taux, tauy, tauz];

    % vectors
    q = [q0; q1; q2; q3];
    w = [wx; wy; wz];
    I_vec = [Ixx, Iyy, Izz];
    B1 = [Bx1; By1; Bz1];
    B2 = [Bx2; By2; Bz2];
    theta = [dt1; dt2; dt3];

    % dynamics
    R_ECI2Body = RotationConversion_sym_v1('QtoDCM', transpose(q));
%     R_ECI2Body = Rq(transpose(q));
    DynOpt.f = attitude_Kin_eqs(transpose(w), q, I_vec, params, R_ECI2Body); 

    %%%%%% global vars store %%%%%%%%
    DynOpt.q_sym = q;
    DynOpt.w_sym = w;

    % map generation
    DynOpt.X = [DynOpt.q_sym; DynOpt.w_sym];
      
    %%% second magnetometer %%%
    DynOpt.Rtheta = RotationConversion_sym_v1('EA321toDCM', transpose(360*theta/pi));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DynOpt.h = DynOpt.w_sym;
    if DynOpt.nMagneto >= 1
        DynOpt.h = [DynOpt.h; R_ECI2Body*B1];
    end
    if DynOpt.nMagneto >= 2
        DynOpt.h = [DynOpt.h; DynOpt.Rtheta*R_ECI2Body*B2];
    end 
     
    %%% linearisation %%%
    % linearization of satellite equations
    DynOpt.G = jacobian(DynOpt.f,DynOpt.X);    
    DynOpt.H = jacobian(DynOpt.h,DynOpt.X); 
    
    DynOpt.Gsym = symfun(DynOpt.G,symarray_G);
    DynOpt.Hsym = symfun(DynOpt.H,symarray_H);
    
    % output
    DynOpt_out = DynOpt;
    params_out = params;

end

