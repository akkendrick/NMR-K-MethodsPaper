function plotWellProfiles_wSonicT2B(figname,waterTable,z,T2dist,T2logbins,...
    DPPdat, SDR_K_all, SDR_K_best, plotNum, fig, sonicCoreT2BData) 
   
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
    title(figname)
    set(gca, 'YDir','reverse')
    %ylim([minDepth,maxDepth])
    ylim([minDepth,maxDepth])

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
    
    plot(DPP_K, DPP_depths,'or','MarkerSize',8)

    %Plot SDR model K estimates
    plot(SDR_K_best, DPP_depths,'.','MarkerSize',24,'Color', [0, 0.4470, 0.7410])

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
    
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])

    set(gca,'XScale','log')

    set(gca,'FontSize',10)

    xlim([10^-6, 10^-2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Plot T2B values
    
    subplot(2,6,plotNum + 2);
    
    box on
    grid on

    hold on
        
    T2BDepths = sonicCoreT2BData.Depthm;
    T2BVals = sonicCoreT2BData.T2Bpeak;
    
    if plotNum <= 6
        plot(T2BVals, T2BDepths,'.','MarkerSize',15,'Color',[0.4660 0.6740 0.1880])
    else
        T2BEst = 2.0680 * 10^3; % Estimate T2B for Wisc based on USGS Data
        % this is T2B from KGM model equation for T2B at 7 deg C estimated from
        % map of groundwater temperature
        % (https://pubs.usgs.gov/wsp/0520f/report.pdf) 
        % "Temperature of water available for industrial use in the United States: 
        % Chapter F in Contributions to the hydrology of the United States,
        % 1923-1924" by W.D. Collins
        % Equation for T2B from Dlugosch et al. 2013
        
        T2BVals = repmat(T2BEst,length(DPP_depths));
        plot(T2BVals, DPP_depths,'.','MarkerSize',15,'Color',[0.4660 0.6740 0.1880])

    end
    
    xlabel(strcat('T_{2B} Estimate', '\rm (ms)'))
    xlim([0,2500])
    
    set(gca, 'YDir','reverse')

    ylim([minDepth,maxDepth])
    ylabel('Depth (m)')
    xlim
    
    set(gca,'FontSize',10)

   
    
end




