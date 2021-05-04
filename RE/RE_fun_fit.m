%% function for the fitting
function out = RE_fun_fit(dn,n,v,d,a,b,version)

    if strcmp(version,'dn')
        out = a*n+b*(v/(2*pi*0.96))*(3e19/d);
    elseif strcmp(version,'n')
        out = a*dn-b/a*(v/(2*pi*0.96))*(3e19/d);
    else
        out = -1;
    end
end