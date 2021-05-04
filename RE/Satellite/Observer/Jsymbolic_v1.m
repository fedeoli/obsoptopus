%% Finding the Magnetic Bias from magnetic measurements

clear all 
close all


nagent = 2
edges  = (nagent-1)*nagent/2
Wuwb = sym('Wuwb','real');
Wgps = sym('Wgps','real');
Wsigma = sym('Wsigma','real');

Chi = sym('Chi%d%d',[nagent,3]);%eval(['Chi' num2str(i) '' num2str(j) ' = sym(''Chi' num2str(i) '' num2str(j) ''');']);
GPS = sym('GPS%d%d',[nagent,3]);%eval(['Chi' num2str(i) '' num2str(j) ' = sym(''Chi' num2str(i) '' num2str(j) ''');']);
Chi0 = sym('Chi0%d%d',[nagent,3]);%eval(['Chi' num2str(i) '' num2str(j) ' = sym(''Chi' num2str(i) '' num2str(j) ''');']);

%UWB distance readings
for i = 1: nagent,
    for j = i+1: nagent,
        eval(['d_' num2str(i) '_' num2str(j) ' = sym(''d_' num2str(i) '_' num2str(j) ''');']);
    end
end

%estimated distances
Juwb = 0;
for i = 1: nagent,
    for j = i+1: nagent,
        eval(['dhat_' num2str(i) '_' num2str(j) ' = sqrt(sum((Chi(' num2str(i) ',:) - Chi(' num2str(j) ',:) ).^2));']);
        eval(['zu_' num2str(i) '_' num2str(j) ' = sym(''zu_' num2str(i) '_' num2str(j) ''');']);
        eval(['Juwb = Juwb + Wuwb*zu_' num2str(i) '_' num2str(j) ' * (dhat_' num2str(i) '_' num2str(j) ' - d_' num2str(i) '_' num2str(j) ')^2;'   ]);
    end
end

Jgps = 0;
for i = 1: nagent,
    eval(['zg_' num2str(i) ' = sym(''zg_' num2str(i)  ''')']);
    eval(['Jgps = Jgps + Wgps*zg_' num2str(i) ' *  (sqrt(sum((Chi(' num2str(i) ',:) - GPS(' num2str(i) ',:) ).^2)))^2'   ]);
end

%initial conditions
Jsigma = 0;
for i = 1: nagent,
    eval(['Jsigma = Jsigma + Wsigma *  (sqrt(sum((Chi(' num2str(i) ',:) - Chi0(' num2str(i) ',:) ).^2)))^2'   ]);
end

J = Juwb + Jgps + Jsigma ;

%Jacobian of J
for i = 1: nagent,
    eval(['gradJ(' num2str(1+3*(i-1)) ') = diff(J,Chi(' num2str(i) ',1)); ']);
    eval(['gradJ(' num2str(2+3*(i-1)) ') = diff(J,Chi(' num2str(i) ',2)) ;']);
    eval(['gradJ(' num2str(3+3*(i-1)) ') = diff(J,Chi(' num2str(i) ',3)); ']);
    
end

%Hessian of J
for i = 1: nagent,
    for j = 1: nagent,
        eval(['HJ(' num2str(1+3*(i-1)) ',' num2str(1+3*(j-1)) ') = diff(gradJ(' num2str(1+3*(i-1)) ') ,Chi(' num2str(j) ',1) );  ']);
        eval(['HJ(' num2str(2+3*(i-1)) ',' num2str(2+3*(j-1)) ') = diff(gradJ(' num2str(2+3*(i-1)) ') ,Chi(' num2str(j) ',2) );  ']);
        eval(['HJ(' num2str(3+3*(i-1)) ',' num2str(3+3*(j-1)) ') = diff(gradJ(' num2str(3+3*(i-1)) ') ,Chi(' num2str(j) ',3) );  ']);
    end
%     for j = 1: nagent,
%         eval(['HJ(' num2str(2+3*(i-1)) ',' num2str(2+3*(j-1)) ') = diff(gradJ(' num2str(2+3*(i-1)) ') ,Chi(' num2str(i) ',' num2str(j) ') ) ; ']);
%     end
%     for j = 1: nagent,
%         eval(['HJ(' num2str(3+3*(i-1)) ',' num2str(3+3*(j-1)) ') = diff(gradJ(' num2str(3+3*(i-1)) ') ,Chi(' num2str(i) ',' num2str(j) ') );  ']);
%     end
end


% The expression are now printed so that can be copy paste into a function that 
% defines the single stem NR-algorithm


%gradJ = simplify(gradJ)
%HJ = simplify(HJ)
%sol = solve(gradJ,Chi)

%
%f = reshape(double(eval(Chi)),3*nagent,1) - pinv(double(eval(HJ)))*(double(eval(gradJ))');
%matlabFunction(gradJ,'File','ObsGradJ')
%matlabFunction(HJ,'File','ObsHJ')





%% % TESTING: PROVIDING A NUMERICAL EXAMPLE
Wuwb = 0;
Wgps = 1;
Wsigma = 0;


for i = 1: nagent,
    for j = 1: 3,
        eval(['X(' num2str(i) ',' num2str(j) ') = ' num2str(randn(1)*i*j)  ';']);
        eval(['GPS' num2str(i) '' num2str(j) ' = X(' num2str(i) ',' num2str(j) ') + ' num2str(1*randn(1))   ';']);
        eval(['Chi0' num2str(i) '' num2str(j) ' = GPS' num2str(i) '' num2str(j)    ';']);
        eval(['Chi' num2str(i) '' num2str(j) ' = Chi0' num2str(i) '' num2str(j) ';']);
    end
end



for i = 1: nagent,
    for j = i+1: nagent,
        eval(['d_' num2str(i) '_' num2str(j) ' = sqrt(sum((X(' num2str(i) ',:) - X(' num2str(j) ',:) ).^2));']);
        eval(['zu_' num2str(i) '_' num2str(j) ' = 0']);
    end
    eval(['zg_' num2str(i) ' = 1']);
end

%%
x(k+1) = x(k) - pinv(Hf)*Jf
MySol = (reshape(double(eval(Chi)),3*nagent,1) - pinv(double(eval(HJ)))*(double(eval(gradJ))') )'
reshape(X,3*nagent,1)'
reshape(X,3*nagent,1)'-MySol

%%
        
        




