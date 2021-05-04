%% symbolic analysys
function [DynOpt_out,params_out] = SymAnalysis_RL_v1

    %%% global vars %%%
    global DynOpt params
    
    %%%%%%%%%%%%%%%%%%%% GENERATE MAPS %%%%%%%%%%%%%%%%%%
    %%%%% Sym analysis %%%%%
    fprintf('Setting Observer parameters\n');
    syms wx wy wz;                      % System angular velocity
    syms q1 q2 q3 q0;                   % Quaternion
    syms Bx1 By1 Bz1 Bx2 By2 Bz2;       % Magnetometer measures

    %%%%% control analysys %%%%%
    syms qrx qry qrz;       % desired attitude

    % vectors
    q = [q0; q1; q2; q3];
    qr = [qrx; qry; qrz];
    w = [wx; wy; wz];

    % numerical data
    I_vec = [params.sat(1).I(1,1); params.sat(1).I(2,2); params.sat(1).I(3,3)];
    I = diag(I_vec);
    r = params.SatellitesCoordinates(1:3);
    r_norm = norm(r);
    v = params.SatellitesCoordinates(4:6);
    v_norm = norm(v);
    mie = params.mi;

    %%%%%%%%%%% quaternion dynamics %%%%%%%%%%%
    % cross product
    Om = [0, -w(1), -w(2), -w(3);
        w(1), 0, w(3), -w(2);
        w(2), -w(3), 0, w(1);
        w(3), w(2), -w(1), 0];
    % f map
    Fq = 0.5*Om*q;

    %%%%%%%%%%% Inertial dynamics %%%%%%%%%%%
%     R_ECI2Body = RotationConversion_sym_v1('QtoDCM', transpose(q));
    R_ECI2Body = Rq(transpose(q));
    Iom = I*w;
    cross_omIom = cross(w,Iom);

    % gravity gradient
    vers_o_Body = -R_ECI2Body*(r/r_norm);
    Io = I*vers_o_Body;
    cross_oIo = cross(vers_o_Body,Io);
    GG_torque = (3*mie/(r_norm^3))*cross_oIo;

    % final dynamics
    Fw = I\( - cross_omIom + GG_torque);
    
    %%%%%% WRAPPED VERSION %%%%
%     F = attitude_Kin_eqs_RL_v1(w, q, I_vec, params, R_ECI2Body); 
%     Fq = F(1:4);
%     Fw = F(5:7);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%% output mapping %%%%%%%%%%%%%
%     R_ECI2Body = RotationConversion_sym_v1('QtoDCM', transpose(q));
    h_temp = w;
    % Magneto measures - Inertial frame
    if DynOpt.nMagneto >= 1
        Magneto1 = [Bx1; By1; Bz1];
        BI_1 = Magneto1;
        Bbody_1 = R_ECI2Body*BI_1;
        h_temp = [h_temp; Bbody_1];
    end
    if DynOpt.nMagneto >= 2
        Magneto2 = [Bx2; By2; Bz2];
        BI_2 = Magneto2;
        Bbody_2 = R_ECI2Body*BI_2;
        h_temp = [h_temp; Bbody_2];
    end

    %%%%%% global vars store %%%%%%%%
    DynOpt.q_sym = q;
    DynOpt.w_sym = w;
    DynOpt.mie_sym = mie;
    DynOpt.v_sym = v;
    DynOpt.v_norm_sym = v_norm;
    DynOpt.qr = qr;
    
    if DynOpt.nMagneto >= 1
        DynOpt.magneto1 = Magneto1;
    end
    if DynOpt.nMagneto >= 2
        DynOpt.magneto2 = Magneto2;
    end
    
    % map generation
    DynOpt.X = [DynOpt.q_sym; DynOpt.w_sym];
    DynOpt.f = [Fq; Fw];
    DynOpt.h = h_temp;

    % output
    DynOpt_out = DynOpt;
    params_out = params;

end

