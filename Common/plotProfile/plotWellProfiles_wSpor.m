function plotWellProfiles_wSpor(waterTable,z,T2dist,T2logbins,...
    DPPdat, SDR_K_all, SDR_K_best, plotNum, fig, coreSampleLoc, sporDepths, sporVal) 
   
    Kline = 5*10^-5;
    

    %T2dist = flip(T2dist);
    minDepth = 1;
    maxDepth = 18;
    
    depths = T2dist(:,1);
    T2dist = T2dist(:,2:end);
   
    figure(fig)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot T2 distribution
    
    % Compute T2ML for all points
    for kk = 1:length(depths)
        T2ML_profile(kk) = sum(T2dist(kk,:).*T2logbins)./sum(T2dist(kk,:));
    end
    
    T2ML_lin = 10.^T2ML_profile;
    T2linbins = 10.^T2logbins;
    
    for kk = 1:length(depths)
            [val, index] = min(abs(T2ML_lin(kk)-T2logbins));
            %Amp_profile(kk) = T2dist(kk,index);
            Amp_profile(kk) = 3; % fixing high value to plot on top
    end
    
    AmpTest = ones(size(z)).*4;
    
    subplot(2,6,plotNum)
    %imagesc(T2dist)
    box on
    grid on
    hold on
    
    waterTableLine = ones(1,length(T2logbins));
    waterTableLine = waterTableLine .* waterTable;
%     for kk = 1:length(depths)
%         currentDepth = depths(kk,:) .* ones(1,100);
%         %plot3(T2logbins,currentDepth,T2dist(kk,:),'b')
%         % ribbon(T2dist(kk,:),currentDepth)
% 
%     end
    
    surf(T2logbins,depths,T2dist,'EdgeColor','none')
    view(0,90)
    %h = pcolor(T2logbins,depths,T2dist);
    %set(h, 'EdgeColor', 'none');
    plot3(T2logbins,waterTableLine,ones(1,length(T2logbins)),'m','LineWidth',3)
    %plot3(log10(T2ML), z, ones(1,length(T2ML)),'*r','MarkerSize',8)
    plot3(T2ML_profile,depths,Amp_profile,'*r','MarkerSize',8)
    
  %  T2ML_test = log10(T2ML);
  %  plot3(T2ML_test,z,AmpTest,'^g','MarkerSize',8)

%     for kk = 1:length(depths)
%         depths(kk)
%         currentDepth = ones(length(T2dist),1)*depths(kk);
%         plot3(T2logbins',currentDepth,T2dist)
%     
%     end
        
    set(gca, 'YDir','reverse')
    %ylim([minDepth,maxDepth])
    ylim([minDepth,maxDepth])
    ylabel('Depth (m)')

    xlim([min(T2logbins),max(T2logbins)])
    xlabel('NMR log_{10} T_2 (s)')

    set(gca,'FontSize',10)
  %   set(gca,'xscale','log')
    
  %  set(gca, 'YDir','reverse')

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Plot K Estimates
    subplot(2,6,plotNum +1);
    
    
    DPP_depths = DPPdat(:,1);
    DPP_K = DPPdat(:,2)*1.16e-5;
    
    hold on
    grid on 
    box on
    
    %scatter(K,z,60,'Filled')

%     [nrow, ncol] = size(k_estimates);
%     for a = 1:ncol
%         plot(k_estimates(:,a),z,k_sym{a},'MarkerSize',5)
%     end
    
    %scatter(dlubacModel,z,40,[0.9290,0.6940,0.1250],'Filled')
    %scatter(SDRModel,z,40,[0.466,0.6740,0.1880],'Filled')
    
    %highDPP = DPP_K(DPP_K > Kline);
    %lowDPP = DPP_K(DPP_K <= Kline);
    
    %highDPP_depths = DPP_depths(DPP_K > Kline);
    %lowDPP_depths = DPP_depths(DPP_K <= Kline);
    
    plot(DPP_K, DPP_depths,'or','MarkerSize',8)

    %plot(highDPP, highDPP_depths,'or','parent',ax1,'MarkerSize',8)
    %plot(lowDPP, lowDPP_depths,'ok','parent',ax1,'MarkerSize',8)
    
    %Plot SDR model K estimates
    plot(SDR_K_best, DPP_depths,'.','MarkerSize',20,'Color', [0, 0.4470, 0.7410])

    
    for j = 1:length(SDR_K_all)
            plot(SDR_K_all{j}, DPP_depths,'.','MarkerSize',4,'Color',[17 17 17]/255)
    end
    
    if plotNum <= 4 
        %plot(6*10^-4, coreSampleLoc,'k*','parent',ax1,'MarkerSize',8)
    end
    
    plot([Kline, Kline], [0,18],'LineWidth',2,'Color',[0.6350, 0.0780, 0.1840])
    box on
    grid on
    
   
    xlabel(strcat('DPP \it K', '\rm (m/s)'))


    set(gca,'XScale','log')

    set(gca,'FontSize',10)
    
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
    ylabel('Depth (m)')

    xlim([10^-6, 10^-2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Plot rho values
    
    subplot(2,6,plotNum + 2);
    
    box on
    grid on

    hold on
        
    if plotNum <= 6
        plot(sporVal,sporDepths,'.','MarkerSize',15)
       
    end
    
    xlabel(strcat('Spor Values', '\rm (1/um)'))
    set(gca, 'YDir','reverse')

    ylim([minDepth,maxDepth])
    ylabel('Depth (m)')
    
    set(gca,'FontSize',10)

    
end




