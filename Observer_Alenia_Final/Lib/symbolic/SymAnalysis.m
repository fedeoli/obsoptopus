%% symbolic analysys
%%%%% global vars %%%%%%
global DynOpt params

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

% differentiation
for i=1:length(Fq)
    for j=1:length(q)
        Fq_q(i,j) = diff(Fq(i),q(j));
    end
end
for i=1:length(Fq)
    for j=1:length(w)
        Fq_w(i,j) = diff(Fq(i),w(j));
    end
end
for i=1:length(Fq)
    for j=1:length(I)
        Fq_theta(i,j) = diff(Fq(i),I(j,j));
    end
end

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
% Fw = I\( - cross_omIom + GG_torque + tau);
Fw = I\( - cross_omIom + GG_torque);

% differentiation
for i=1:length(Fw)
    for j=1:length(q)
        Fw_q(i,j) = diff(Fw(i),q(j));
    end
end
for i=1:length(Fw)
    for j=1:length(w)
        Fw_w(i,j) = diff(Fw(i),w(j));
    end
end
for i=1:length(Fw)
    for j=1:length(I)
        Fw_theta(i,j) = diff(Fw(i),I(j,j));
    end
end

%%%%%%%%%%%%%%%% final store %%%%%%%%%%%%%
dF = [Fq_q, Fq_w, Fq_theta;
      Fw_q, Fw_w, Fw_theta];
  
dF = simplify(dF,'IgnoreAnalyticConstraints',true);
dF = simplifyFraction(dF);
  
%%%%%% global vars store %%%%%%%%
DynOpt.q_sym = q;
DynOpt.w_sym = w;
DynOpt.I_sym = I;
DynOpt.I_vec_sym = I_vec;
DynOpt.r_sym = r;
DynOpt.r_norm_sym = r_norm;
DynOpt.mie_sym = mie;
DynOpt.v_sym = v;
DynOpt.v_norm_sym = v_norm;
DynOpt.qr = qr;

DynOpt.dF_sym = dF;
DynOpt.Fq_sym = Fq;
DynOpt.Fw_sym = Fw;

