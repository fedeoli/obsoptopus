%{
How to use the program:
Suppose the function to be evluated is:
f(x)=6+3cos(x)
with an interval of 0 and pi/2

Enter the function to be integrated: 6+(3*cos(x))
Enter the number of segments(2 is the minimum):7

*This is arbitrary. 2 segment is a single application of 1/3 Simpson's rule
,3 segment is a single application of 3/8 Simpson's rule. An even # of
segments greater than 3 is a multiple application of 1/3 Simpson's rule. An
odd # of segments greater than 3 is a combination of 1/3 and 3/8 Simson's 
rules.*

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
n=input("Enter the number of segments(2 is the minimum):");
if n<2
    n=0;
    fprintf("Wrong input! The number of segments must be atleast 2. Please try again.\n");
end
end

a=input('The lower limit is:');
b=input('The upper limit is:');
    
%This is for initialization
fa=0;
f4aprime=0;
h=(b-a)/n;
m1=a;
 
%Evaluate the 2nd derivative of the function for the approximate absolute
%error
syms x
d4y=diff(diff(diff(diff(f(x),x))))
f4prime=(int(d4y,a,b))/(b-a);

if (mod(n,2)==0)
    fprintf('Since there are even number of segments,\nthe function is integrated using the 1/3 Simspons rule:\n\n');
    for i=1:n/2
        m2=m1+h;
        m3=m2+h;
        I=f(m1)+(4*f(m2))+f(m3);
        I=I+fa;
        fa=I;
        m1=m3;
    end
    I=(h/3)*I;
    e=-(1/(180*n^4))*(f4prime)*(b-a)^5;
    e=vpa(subs(e));
    fprintf('The solved integral is %.6f with an approximate absolute error of %.6f\n',I,e);
end

if n==3
    fprintf('Since there are three segments,\nthe function is integrated using the 3/8 Simspons rule:\n\n');
    m1=a+h;
    m2=m1+h;
    I=(3*h/8)*(f(a)+(3*f(m1))+(3*f(m2))+f(b));
    e=(-1/6480)*(f4prime)*(b-a)^5;
    e=vpa(subs(e));
    fprintf('The solved integral is %.6f with an approximate absolute error of %.6f\n',I,e);
end


if (n>3)&&(mod(n,2)~=0)
    fprintf('Since the segments is >3 and is odd,\nthe function is integrated using a combination of 1/3 and 3/8 Simspons rule:\n\n');
    n1=n-3;
    for i=1:n1/2
        m2=m1+h;
        m3=m2+h;
        I1=f(m1)+(4*f(m2))+f(m3);
        I1=I1+fa;
        fa=I1;
        m1=m3;
    end
    I1=(h/3)*I1;
    f4prime=(int(d4y,a,m1))/(m1-a);
    e1=-(1/(180*n1^4))*(f4prime)*(m1-a)^5;
    e1=vpa(subs(e1));
   
    m2=m1+h;
    m3=m2+h;
    I2=(3*h/8)*(f(m1)+(3*f(m2))+(3*f(m3))+f(b));
    f4prime=(int(d4y,m1,b))/(b-m1);
    e2=(-1/6480)*(f4prime)*(b-m1)^5;
    e2=vpa(subs(e2));
    
    I=I1+I2;
    e=e1+e2;
    fprintf('The solved integral is %.6f with an approximate absolute error of %.6f\n',I,e);
end