
Nfigure1 = 40;
Nfigure2 = 41;
MyFontSize = 18
kk = 0;
minmin1 = 100E3;
minmin1_kk = 0;
minmin2 = 100E3;
minmin2_kk = 0;

if(ObserverTest.TotalRuns == 1)
    
    figure(Nfigure2)
    errorbar(1:1:ObserverTest.Nagents, ObserverTest.MeanDistances(1,:)*1000,ObserverTest.SigmaDistances(1,:)*1000,'b--','LineWidth',2);
    ylabel('Relative Distances: mean+-sigma [m]','FontSize',MyFontSize,'Interpreter','Latex')
    xlabel('Agents','FontSize',MyFontSize,'Interpreter','Latex')
    set(gca,'FontSize',MyFontSize);
    grid on
    
    figure(Nfigure1)
    errorbar(1:1:ObserverTest.Nagents, ObserverTest.Mean(1,:)*1000,ObserverTest.Sigma(1,:)*1000,'LineWidth',2);
    ylabel('Position: mean+-sigma [m]','FontSize',MyFontSize,'Interpreter','Latex')
    xlabel('Agents','FontSize',MyFontSize,'Interpreter','Latex')
    set(gca,'FontSize',MyFontSize);
    grid on
    
else
    
    
    for np=1:ObserverTest.N_P,
        for nr=1:ObserverTest.N_R,
            
            for j=1:length(ObserverTest.N_Wuwb_range),
                Mean_surf1 = zeros(length(ObserverTest.N_Wgps_range),length(ObserverTest.N_Wsigma_range));
                Mean_surf2 = zeros(length(ObserverTest.N_Wgps_range),length(ObserverTest.N_Wsigma_range));
                for ngps=1:length(ObserverTest.N_Wgps_range),
                    for nsigma=1:length(ObserverTest.N_Wsigma_range),
                        kk = kk+1;
                        if(min(ObserverTest.UnfeasableR(kk,:))<0.1)
                            Mean_surf1(ngps,nsigma) = mean(ObserverTest.Mean(kk,:))*1E3;
                            Mean_surf2(ngps,nsigma) = mean(ObserverTest.MeanDistances(kk,:))*1E3;
                            if(minmin1 > Mean_surf1(ngps,nsigma)+mean(ObserverTest.Sigma(kk,:)) )
                                minmin1 = Mean_surf1(ngps,nsigma) + mean(ObserverTest.Sigma(kk,:));
                                minmin1_kk = kk;
                                ObserverTest.Gains(:,minmin1_kk)
                            end
                            if(minmin2 > Mean_surf2(ngps,nsigma)+mean(ObserverTest.SigmaDistances(kk,:)) )
                                minmin2 = Mean_surf2(ngps,nsigma) + mean(ObserverTest.SigmaDistances(kk,:));
                                minmin2_kk = kk;
                                ObserverTest.Gains(:,minmin2_kk)
                            end
                        else
                            Mean_surf1(ngps,nsigma) = 0;
                            Mean_surf2(ngps,nsigma) = 0;
                        end
                    end
                end
                grid on
                if(length(ObserverTest.N_Wgps_range)==1)
                    figure(Nfigure1)
                    plot(ObserverTest.N_Wsigma_range,Mean_surf1,'LineWidth',2);
                    ylabel(['[m] W_{uwb} =  ' num2str(ObserverTest.N_Wuwb_range(j)) ',  Ri_{1,1}=' num2str(eig(ObserverTest.RiAll(kk,1))) ', Pi_{1,1} = ' num2str(eig(ObserverTest.PiAll(kk,1))) ],'FontSize',MyFontSize,'Interpreter','Latex')
                    grid on
                    figure(Nfigure2)
                    plot(ObserverTest.N_Wsigma_range,Mean_surf2,'LineWidth',2);
                    ylabel([' [m] W_{uwb} =  ' num2str(ObserverTest.N_Wuwb_range(j)) ',  Ri_{1,1}=' num2str(eig(ObserverTest.RiAll(kk,1))) ', Pi_{1,1} = ' num2str(eig(ObserverTest.PiAll(kk,1))) ],'FontSize',MyFontSize,'Interpreter','Latex')
                    grid on
                else
                    figure(Nfigure1)
                    surf(ObserverTest.N_Wsigma_range,ObserverTest.N_Wgps_range,Mean_surf1);
                    zlabel(['[m] W_{uwb} =  ' num2str(ObserverTest.N_Wuwb_range(j)) ',  Ri_{1,1}=' num2str(eig(ObserverTest.RiAll(kk,1))) ', Pi_{1,1} = ' num2str(eig(ObserverTest.PiAll(kk,1))) ],'FontSize',MyFontSize,'Interpreter','Latex')
                    ylabel(['W_{gps}'],'FontSize',MyFontSize,'Interpreter','Latex');
                    grid on
                    figure(Nfigure2)
                    surf(ObserverTest.N_Wsigma_range,ObserverTest.N_Wgps_range,Mean_surf2);
                    zlabel(['[m] W_{uwb} =  ' num2str(ObserverTest.N_Wuwb_range(j)) ',  Ri_{1,1}=' num2str(eig(ObserverTest.RiAll(kk,1))) ', Pi_{1,1} = ' num2str(eig(ObserverTest.PiAll(kk,1))) ],'FontSize',MyFontSize,'Interpreter','Latex')
                    ylabel(['W_{gps}'],'FontSize',MyFontSize,'Interpreter','Latex');
                    grid on
                end
                xlabel(['W_{sigma}'],'FontSize',MyFontSize,'Interpreter','Latex');
                if(ObserverTest.CentralizedOptimization==1)
                    if(ObserverTest.UWB_StepsBackInTime == 0)
                        figure(Nfigure1)
                        title(['Position - Centralized W_{GPSall}=' num2str(ObserverTest.Weight_GPSall)],'FontSize',MyFontSize,'Interpreter','Latex');
                        figure(Nfigure2)
                        title(['Distances - Centralized W_{GPSall}=' num2str(ObserverTest.Weight_GPSall)],'FontSize',MyFontSize,'Interpreter','Latex');
                    else
                        figure(Nfigure1)
                        title(['Window - Position - Centralized W_{GPSall}=' num2str(ObserverTest.Weight_GPSall)],'FontSize',MyFontSize,'Interpreter','Latex');
                        figure(Nfigure2)
                        title(['Window - Distances - Centralized W_{GPSall}=' num2str(ObserverTest.Weight_GPSall)],'FontSize',MyFontSize,'Interpreter','Latex');
                    end
                else
                    if(ObserverTest.UWB_StepsBackInTime == 0)
                        figure(Nfigure1)
                        title(['Position - Decentralized'],'FontSize',MyFontSize,'Interpreter','Latex');
                        figure(Nfigure2)
                        title(['Distances - Decentralized'],'FontSize',MyFontSize,'Interpreter','Latex');
                    else
                        figure(Nfigure1)
                        title(['Window - Position - Decentralized'],'FontSize',MyFontSize,'Interpreter','Latex');
                        figure(Nfigure2)
                        title(['Window - Distances - Decentralized'],'FontSize',MyFontSize,'Interpreter','Latex');
                    end
                end
                %if(np<= ObserverTest.N_P && nr<ObserverTest.N_R)
                pause
                %end
            end
        end
    end
end
