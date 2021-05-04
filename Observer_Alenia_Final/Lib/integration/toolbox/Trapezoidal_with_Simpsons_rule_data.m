%{
How to use the program:
Suppose the data to be evluated is:
x:    0.5	  2	      3	      4       6       8      10	     11
f(x): 336	294.4	266.4	260.8	260.5	249.6	193.6	165.6

Input the x data in one-column matrix form used in matlab:
[0.5 2 3 4 6 8 10 11]'

Input the y data in one-column matrix form used in matlab:
[336 294.4 266.4 260.8 260.5 249.6 193.6 165.6]'
%}

%This is for inputs
X=input('Input the x data in one-column matrix form used in matlab:')
Y=input('Input the y data in one-column matrix form used in matlab:')           
n=size(X,1);
m=size(Y,1);
fa=0;
s=0;
a=0;
d=0;

%This is the main loop program

if (n==3)
    h1=X(2)-X(1);
    h2=X(3)-X(2);
    if abs(h1-h2)<0.00001
       I=(h1/3)*(Y(1)+(4*Y(2))+Y(3));
    end
    
    if abs(h1-h2)>0.00001
       for i=1:n-1
           h=X(i+1)-X(i);
           w=(Y(i+1)+Y(i))/2;
           I=fa+(h*w);
           fa=I;
       end
    end
end

if (n==2)
    I=(X(2)-X(1))*(Y(2)+Y(1))/2;
end

for i=4:n
    if (s>a)
        s=s-1;
        continue; 
    end
    
    a=s;
    h1=X(i-2)-X(i-3);
    h2=X(i-1)-X(i-2);
    h3=X(i)-X(i-1);
   
    if (abs(h1-h2)>0.00001)
        I=fa+(h1*(Y(i-3)+Y(i-2))/2);
    end   
    
    if (abs(h1-h2)<0.00001)&&(abs(h2-h3)>0.00001)
       I=fa+((2*h1/6)*(Y(i-3)+(4*Y(i-2))+Y(i-1)));
       s=s+1;  
    end
    
    if (abs(h1-h2)<0.00001)&&(abs(h2-h3)<0.00001)
       I=fa+((3*h1/8)*(Y(i-3)+(3*Y(i-2))+(3*Y(i-1))+Y(i)));
       s=s+2;  
       if i==n-1
           d=1;
       end
    end
    fa=I;

end
if (a==s)
    h1=X(i-1)-X(i-2);
    h2=X(i)-X(i-1);
    
    if (abs(h1-h2)>0.00001)
    I=fa+(h1*(Y(i-1)+Y(i-2))/2);
    I=I+(h2*(Y(i)+Y(i-1))/2);  
    end
    
    if (abs(h1-h2)<0.00001)
    I=fa+((2*h1/6)*(Y(i-2)+(4*Y(i-1))+Y(i)));    
    end
end

if (a==s-1)
    h3=X(i)-X(i-1);
    I=fa+(h3*(Y(i)+Y(i-1))/2); 
end
%Interpretation of results
fprintf('\nThe solved integral is %.6f\n',I);