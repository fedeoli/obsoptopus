%% plot RL results
function plot_RL(RL)

        %%% RL results 
        figure
        b = bar3(RL.S.Q');
        set(gca,'XTickLabel',RL.A.domain_amp_grid)
        xlabel('Amp')
        set(gca,'YTickLabel',RL.A.domain_d_grid)
        ylabel('Duty')
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
        end
        
        figure
        b = bar3(RL.S.N);
        set(gca,'XTickLabel',RL.A.domain_amp_grid)
        xlabel('Amp')
        set(gca,'YTickLabel',RL.A.domain_d_grid)
        ylabel('Duty')
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
        end
        
        %%% PWM results
        figure
        for i = 1:length(RL.S.pos_best)
            hold on
            plot(RL.S.DynOpt(RL.S.pos_best(i)).True_quat(1,:))
        end
        plot(RL.S.DynOpt(RL.S.pos_best(1)).desatt_true(1,:),'--')
        
        figure
        for i = 1:length(RL.S.pos_compare)
            hold on
            plot(RL.S.DynOpt(RL.S.pos_compare(i)).True_quat(1,:))
        end
        plot(RL.S.DynOpt(RL.S.pos_compare(1)).desatt_true(1,:),'--')
        
        %%% eigenvalues results
%         figure
%         for i = 1:length(RL.S.pos_best)
%             hold on
%             plot(RL.S.DynOpt(RL.S.pos_best(i)).dtheta_eig_min_normalised)
%         end
%         
%         figure
%         for i = 1:length(RL.S.pos_compare)
%             hold on
%             plot(RL.S.DynOpt(RL.S.pos_compare(i)).dtheta_eig_min_normalised)
%         end
        
end