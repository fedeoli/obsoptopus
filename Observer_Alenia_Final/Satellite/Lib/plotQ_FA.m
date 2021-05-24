%% function plot Q
function Qval = plotQ_FA(Q,domain_S,domain_A,N_count)

    
    %%% show results
    figure(1)
    clf
    a = max(Q,[],2);
    Qval = [];
    for j=1:size(a,1)
       Aval = find(Q(j,:) == a(j));
       str_array = string(1+domain_A(:,Aval)');
       Qval(j) = strcat(str_array(1), str_array(2));
    end
    
    size_Q = [size(domain_S,2)-1,size(domain_S,2)-1];
    Qval = reshape(Qval,size_Q);
    heatmap(domain_S(2,2:end),domain_S(1,2:end),Qval);
    xlabel('edot')
    ylabel('e')
    
    % legend
    title({'1 = disabled, 2 = enabled','first = magneto, second = input'}) 

    % workaround - x on top
    ax = gca;
    axp = struct(ax);
    axp.Axes.XAxisLocation = 'top';
    
    %%% N count %%%
    figure(2)
    heatmap(domain_S(2,2:end),domain_S(1,2:end),N_count)
    xlabel('edot')
    ylabel('e')
    % workaround - x on top
    ax = gca;
    axp = struct(ax);
    axp.Axes.XAxisLocation = 'top';
    % legend
    title('N visits')  
    
    %%% comnined view %%%
    figure(3)
    heatmap(domain_S(2,2:end),domain_S(1,2:end),sign(N_count).*Qval)
    xlabel('edot')
    ylabel('e')
    % workaround - x on top
    ax = gca;
    axp = struct(ax);
    axp.Axes.XAxisLocation = 'top';
    % legend
    title('Combined view')  
    
end