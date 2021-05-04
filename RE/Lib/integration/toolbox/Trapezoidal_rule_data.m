%{
How to use the program:
Suppose the data to be evluated is:
x:    -2	0	2	 4	 6	 8	 10
f(x): 35	5	-10	 2	 5	 3	 20

Input the x data in one-column matrix form used in matlab:
[-2 0 2 4 6 8 10]'

Input the y data in one-column matrix form used in matlab:
[35 5 -10 2 5 3 20]'
%}

%This is for inputs
X=input('Input the x data in one-column matrix form used in matlab:')
Y=input('Input the y data in one-column matrix form used in matlab:')           
n=size(X);
m=size(Y);
fa=0;

%This is the main loop program
for i=1:n-1
    h=X(i+1)-X(i);
    w=(Y(i+1)+Y(i))/2;
    I=fa+(h*w);
    fa=I;
end


%Interpretation of results
fprintf('\nThe solved integral is %.6f\n',I);