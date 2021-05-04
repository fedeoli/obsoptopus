%{
How to use the program:
Suppose the function to be evaluated is:
f(x)= -0.0547x^4 +0.8646x^3 -4.1562x^2 +6.2917x +2
to be integrated from 0 to 8 with 0.10% tolerance oor upto 3 iterations

Enter the function to be integrated: 
(-0.0547*x^4)+(0.8646*x^3)-(4.1562*x^2)+(6.2917*x)+2

The lower limit is:0
The upper limit is:8

If tolerance of 0.10%:
Please choose what termination criteria you desire.
Press "1" for approximate relative error and "2" for number of iterations:1
The desired approximate relative error in % is:0.1

If 3 iterations:
Please choose what termination criteria you desire.
Press "1" for approximate relative error and "2" for number of iterations:2
The desired number of iterations is:3

*Iteration 1 has a final value with local truncation error O(h^4). Iteration 2
has a final value with local truncation error of 0(h^6); iteration 3 O(h^8), and
so on.*
%}

%This is for inputs
m=inputdlg('Enter the function to be integrated:');   
s=m{:};                                         
d=str2func(['@(x)' s]);                         
f=str2func(['@(x)' vectorize(s)])               
a=input('The lower limit is:');
b=input('The upper limit is:');
term=0;
A=zeros(10,10);

while (term<1)||(term>2)
    term=input('Please choose what termination criteria you desire.\nPress "1" for approximate relative error and "2" for number of iterations:');
    if term==1
        e=input('The desired approximate relative error in % is:');
        itr=NaN;
        break;
    end
    if term==2
        itr=input('The desired number of iterations is:');
        e=NaN;
        break;
    end
    fprintf('Wrong input! Please try again.\n\n');   
end

%This is for initialization






%This is the main loop program
for l=1:100
    fprintf('\nIteration %.0f:\n\n',l);
    r=l;
    n=1;
    I=zeros(l+1,1);
  
    for j=1:l+1
        h=(b-a)/n;
        X=a;
        fa=0;
        for i=1:n+1
            fX=2*f(X);
            fX=fa+fX;
            fa=fX;
            X=X+h;   
        end
        fX=fX-f(a)-f(b);
        I(j,1)=(h*fX)/2;
        A(j,1)=I(j,1);
        n=2^j;
    end
    
    for k=2:l+1
        if k>2
           r=r-1;
        end
            
        for m=1:r
        I(m,1)=(((4^(k-1))*(I(m+1,1)))-(I(m,1)))/((4^(k-1))-1);
        A(m,k)=I(m,1);
        g=abs((I(m,1)-I(m+1,1))/I(m,1))*100;
        end
    end
    
    A=A(1:k,1:k);
    disp(A);
    fprintf('The percent relative error is ');
    disp(g);
    if (l==itr)||(g<e)
        break;
    end
end 

%Interpretation of results
fprintf('After %.0f iterations, the solved integral is %.7f with an approximate\nrelative error of %.7f and truncation error O(h^%.0f)\n',l,A(1,k),g,2*k);
return
