function J = qopt_1(x)

global qopt_p

Rx = quat2dcn(x);

J = norm(qopt_p.MRx)^2

J = 0;

