e1 = [pi/4,-pi/3,0];
de = [pi/6,pi/10,-pi/4]
q1 = eul2quat(e1); %qt
M1 = quat2dcm(q1);
M2 = angle2dcm(de(1),de(2),de(3))*M1;
q2 = dcm2quat(M2);
dq = quatnormalize(quatmultiply(quatinv(q1),q2))
quat2eul(dq)
 %dehat = quat2eul(dq) - de
%%

M1 = quat2dcm(q1);
M2 = quat2dcm(q2);
dM = M2'*M1;
dqm = dcm2quat(dM);
deq = quat2angle(dqm)

%%


de = 0.05;
N = 100;
dEmax = pi/2;
q = zeros(N,4); 
e = zeros(N,3); 
detM = zeros(N,1); 
esmooth = e; 
for k=1:N,
    eul = [pi/2,-pi/2.1,0+de*k];
    esmooth(k,:) = esmooth(max(1,k-1),:);
    q(k,:) = quatnormalize(eul2quat(eul));
    detM(k) = det(angle2dcm(eul(1),eul(2),eul(3)));
    e(k,:) = quat2eul(q(k,:));
    for j=1:3,
        if( abs(e(k,j) - esmooth(k,j)) > dEmax )
            esmooth(k,j) = e(k,j) - sign(e(k,j) - esmooth(k,j))*2*pi;
        else
            esmooth(k,j) = e(k,j);
        end
    end
end

figure(1)
subplot(3,1,1)
plot(q)
subplot(3,1,2)
plot(e*180/pi)
hold on
plot(esmooth*180/pi,'--')
subplot(3,1,3)
plot(detM)




%%
v = [1;2;-1]
R = quat2dcm(q)
M = R*v
R'*M-v
det(R)
R*R'
clear R
%%
R = sym('R%d%d',[3,3])


%%
syms n x y z b1 b2 b3 dx dy dz
clear n
n = sqrt(1-(x^2+y^2+z^2));
R = [[2*(n^2 + x^2) - 1, 2*(x*y - n*z),  2*(x*z + n*y)]
[2*(x*y + n*z), 2*(n^2 + y^2) - 1, 2*(y*z - n*x)]
[2*(x*z - n*y),  2*(y*z + n*x), 2*(n^2 + z^2) - 1]];

M0 = [1; 2; -1];
M0 = M0/norm(M0);
r1 = pi/6;
r2 = pi/12;
r3 = -pi/3;
q = angle2quat(r1,r2,r3);

qbis = angle2quat(r1-pi/5,r2+pi/5,r3+pi/5);
M = quat2dcm(q)*M0;
Mbis = quat2dcm(qbis)*M0;
n = sqrt(1-((x+dx)^2+(y+dy)^2+(z+dz)^2));
Rbis = [[2*(n^2 + (x+dx)^2) - 1, 2*((x+dx)*(y+dy) - n*(z+dz)),  2*((x+dx)*(z+dz) + n*(y+dy))]
[2*((x+dx)*(y+dy) + n*(z+dz)), 2*(n^2 + (y+dy)^2) - 1, 2*((y+dy)*(z+dz) - n*(x+dx))]
[2*((x+dx)*(z+dz) - n*(y+dy)),  2*((y+dy)*(z+dz)+ n*(x+dx)), 2*(n^2 + (z+dz)^2) - 1]];


f = R*M0;
fbis = (R+ Rbis)*M0;

Jf = jacobian(f,[x,y,z]);
detJf = det(Jf)
sol = solve(f,[x,y,z])

Jfbis  = jacobian(fbis,[x,y,z]);
detJfbis = det(Jfbis)
%solbis = solve(fbis,[x,y,z]) %it  does not find it
optionsObs = optimoptions(@fminunc,  'MaxIter', 100);%'display','off'); % 'Algorithm', 'trust-region','SpecifyObjectiveGradient', true,                   
global qopt_p
qopt_p.M = M;
qopt_p.Mbis = Mbis;
qopt_p.M0 = M0;
sol = fmin(@qopt_1,[0,0,0],optionsObs);


     





%% play with symbolic

clear all
close all

a1 = sym('a1','real')
a2 = sym('a2','real')
p1 = sym('p1','real')
p2 = sym('p2','real')
p3 = sym('p3','real')
c1 = sym('c1','real')
c2 = sym('c2','real')

A = [0 1; -a1 -a2];
P = [p1,p2; p2, p3 ];
temp = A'*P+P*A+[c1,0;0,c2];
%%
sol = solve(temp,[p1,p2,p3])
%% 
p1 = sol.p1
p2 = sol.p2
p3 = sol.p3
P = [p1,p2; p2, p3 ];
S1 = P(1,1)
S2 = det(P)

 
     



