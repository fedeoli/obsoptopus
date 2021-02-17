%{
How to use the program:
Suppose the function to be evluated is:
f(x)=6+3cos(x)
with an interval of 0 and pi/2

Enter the function to be integrated: 6+(3*cos(x))
Enter the number of segments(1 is the minimum):4
*This is arbitrary. 1 segment is a single application of Trapezoidal rule*

The lower limit is:0
The upper limit is:pi/2
%}

%This is for inputs
m=inputdlg('Enter the function to be integrated:');   
s=m{:};                                         
d=str2func(['@(x)' s]);                         
f=str2func(['@(x)' vectorize(s)])               
n=0;
while n==0    
n=input("Enter the number of segments(1 is the minimum):");
if n>=1
    break
end
if n<1
    n=0;
    fprintf("Wrong input! The number of segments must be atleast 2. Please try again.\n");
end
end

a=input('The lower limit is:');
b=input('The upper limit is:');

%This is for initialization
X=a;
h=(b-a)/n;
fa=0;
f2aprime=0;

%Evaluate the 2nd derivative of the function for the approximate absolute
%error
syms x
ddy=diff(diff(f(x),x));

%This is the main loop program
for i=1:n+1
    oldX=X;
    fX=2*f(X);
    fX=fa+fX;
    fa=fX;
    X=X+h;
    f2prime=(int(ddy,oldX,X))/(X-oldX);
    f2prime=f2aprime+f2prime;
    f2aprime=f2prime;
end

fX=fX-f(a)-f(b);
I=(h*fX)/2;
f2prime=f2prime/n;
e=(-1/(12*n^2))*f2prime*(b-a)^3;
e=vpa(subs(e));

%Interpretation of results
fprintf('\nThe solved integral is %.6f with an approximate absolute error of %.4f\n',I,e);