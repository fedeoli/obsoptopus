%% symbolic analysys
function [DynOpt_out,params_out] = SymAnalysis_RL_v2

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

    %%%%% control analysys %%%%%
    syms qrx qry qrz;       % desired attitude

    % vectors
    q = [q0; q1; q2; q3];
    qr = [qrx; qry; qrz];
    w = [wx; wy; wz];
    I = diag([Ixx Iyy Izz]);
    I_vec = [Ixx; Iyy; Izz];
    r = [rx; ry; rz];
    v = [vx; vy; vz];

    %%%%%%%%%%% quaternion dynamics %%%%%%%%%%%
    % cross product
    Om = [0, -w(1), -w(2), -w(3);
        w(1), 0, w(3), -w(2);
        w(2), -w(3), 0, w(1);
        w(3), w(2), -w(1), 0];
    % f map
    Fq = 0.5*Om*q;

    %%%%%%%%%%% Inertial dynamics %%%%%%%%%%%
    R_ECI2Body = RotationConversion_sym_v1('QtoDCM', transpose(q));
    Iom = I*w;
    cross_omIom = cross(w,Iom);

    % gravity gradient
    vers_o_Body = -R_ECI2Body*(r/r_norm);
    Io = I*vers_o_Body;
    cross_oIo = cross(vers_o_Body,Io);
    GG_torque = (3*mie/(r_norm^3))*cross_oIo;

    %%%%%%%%%%%%% control dynamics %%%%%%%%%%%
    R_ECI2Hill = RECI2Hill_sym_v1(r, v);
    omega_Hill2ECI_ECI = OmegaHill2ECI(r, v, mie);
    % Convert the desired attitude to dcm representation
    R_Hill2Body = RotationConversion_sym_v1('EA321toDCM', transpose(qr)*180/pi);

    % Refer the desired attitude wrt ECI reference frame and convert to quaternions
    R_ECI2Body = R_Hill2Body*R_ECI2Hill;

    % convert reference trajectory
    q_ref = RotationConversion_sym_v1('EA321toQ', transpose(qr)*180/pi);

    % testing
    omega_ref = R_ECI2Body*omega_Hill2ECI_ECI;

    % Compute the deputy current angular velocity wrt to Hill reference frame
    omega_Body2Hill_Body = w - omega_ref;

    % Compute the quaternion error between current and desired attitudes
    q_err = QuaternionError_sym_v1(transpose(q), q_ref);

    % Compute PD attitude control torque
    % sign_flag = q_err(1)/abs(q_err(1));
    tau = -sign(q_err(1))*params.sat(1).kp*transpose(q_err(2:4)) - params.sat(1).kd*omega_Body2Hill_Body;


    % final dynamics
    Fw = I\( - cross_omIom + GG_torque + tau);
%     Fw = I\( - cross_omIom + GG_torque);

    % numerical data
    I_vec_num = [params.sat(1).I(1,1); params.sat(1).I(2,2); params.sat(1).I(3,3)];
    I_num = diag(I_vec_num);
    r_num = params.SatellitesCoordinates(1:3);
    r_norm_num = norm(r_num);
    v_num = params.SatellitesCoordinates(4:6);
    v_norm_num = norm(v_num);
    mie_num = params.mi;

    %%%%%% global vars store %%%%%%%%
    DynOpt.q_sym = q;
    DynOpt.w_sym = w;
    DynOpt.mie_sym = mie;
    DynOpt.v_sym = v;
    DynOpt.v_norm_sym = v_norm;
    DynOpt.qr = qr;

    % map generation
    DynOpt.X = [DynOpt.q_sym; DynOpt.w_sym];
    DynOpt.f = [Fq; Fw];
    DynOpt.h = DynOpt.w_sym;

    % output
    DynOpt_out = DynOpt;
    params_out = params;

end

