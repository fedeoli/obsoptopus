function dXdt = mRiccati(t, X, A, Q)
%Convert from "n^2"-by-1 to "n"-by-"n"
X = reshape(X, size(A));

%Determine derivative
dXdt = A.'*X + X*A + Q; 

%Convert from "n"-by-"n" to "n^2"-by-1
dXdt = dXdt(:); 