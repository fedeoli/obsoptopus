function Xnew  = PlantJumpMap(X,j,forwardpropagation,estimate)
% X: state and parameters
% j: index used to retrieve the data (input) on the moving orizon window
% forwardpropagation: 1 if forward propagation, -1 otherwise
% estimate: 1 if used for the estimated plant, 0 othwerwise
% Example: x = A*x+B*u(j); j-1 is the index used to select input value

global DynOpt

x = X(1:DynOpt.StateDim);

% if(estimate==1)
%     A = [0, 1/DynOpt.Ts,0,0; 0,0,1,0; ...
%         0,-DynOpt.c0+(X(6)*max(DynOpt.U(2,j)-1,0)/(DynOpt.Zf^2-X(2)^2)),...
%         -DynOpt.c1,X(5)/(DynOpt.Zh^2-X(2)^2);-DynOpt.Ki*DynOpt.beta,-DynOpt.Kp*DynOpt.beta,-DynOpt.Kd*DynOpt.beta,-DynOpt.beta];
% else
%      A = [0, 1/DynOpt.Ts,0,0; 0,0,1,0; ...
%          0,-DynOpt.c0+(DynOpt.gamma0*max(DynOpt.U(2,j)-1,0)/(DynOpt.Zf^2-X(2)^2)),...
%         -DynOpt.c1,DynOpt.gamma1/(DynOpt.Zh^2-X(2)^2);-DynOpt.Ki*DynOpt.beta,-DynOpt.Kp*DynOpt.beta,-DynOpt.Kd*DynOpt.beta,-DynOpt.beta];
% end
if(estimate==1)
    A = [0, 1,0,0; 0,0,1,0; ...
        0,-DynOpt.c0+X(5),-DynOpt.c1,X(6); %0,-DynOpt.c0+X(6)*max(DynOpt.U(2,j)-1.15,0),-DynOpt.c1,X(5);
        -DynOpt.Ki*DynOpt.gainPHD*DynOpt.beta,-DynOpt.Kp*DynOpt.gainPHD*DynOpt.beta,-DynOpt.Kd*DynOpt.gainPHD*DynOpt.beta,-DynOpt.beta];
else
     A = [0, 1,0,0; 0,0,1,0; ...
         0,-DynOpt.c0+DynOpt.gamma0,-DynOpt.c1,DynOpt.gamma1;%0,-DynOpt.c0+DynOpt.gamma0*max(DynOpt.U(2,j)-1.15,0),-DynOpt.c1,DynOpt.gamma1;
         -DynOpt.Ki*DynOpt.gainPHD*DynOpt.beta,-DynOpt.Kp*DynOpt.gainPHD*DynOpt.beta,-DynOpt.Kd*DynOpt.gainPHD*DynOpt.beta,-DynOpt.beta];
end


if(forwardpropagation == 1)
    %Xnew = x + ones(4,1);%x + DynOpt.Ts*(A*x); 
    x = x + DynOpt.Ts*(A*x); 
else
    %Xnew = x - ones(4,1);%pinv(eye(4)+DynOpt.Ts*A)*x; 
    x = pinv(eye(4)+DynOpt.Ts*A)*x; 
end

Xnew = [x;X(DynOpt.StateDim+1:end)];