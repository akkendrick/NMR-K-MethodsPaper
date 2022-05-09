function plotWellProfiles_stacked(sites,DPPdat)
% Plot 
    
for jj = 1:length(sites)

    %T2dist = flip(T2dist);
    minDepth = 1;
    maxDepth = 18;
    
    

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Plot K Estimates

    DPP_depths = DPPdat{jj}(:,1);
    DPP_K = DPPdat{jj}(:,2)*1.12e-5;
    
    hold on
    
    %scatter(K,z,60,'Filled')

%     [nrow, ncol] = size(k_estimates);
%     for a = 1:ncol
%         plot(k_estimates(:,a),z,k_sym{a},'MarkerSize',5)
%     end
    
    %scatter(dlubacModel,z,40,[0.9290,0.6940,0.1250],'Filled')
    %scatter(SDRModel,z,40,[0.466,0.6740,0.1880],'Filled')
   
    
    plot(DPP_K, DPP_depths,'o','MarkerSize',8,'LineWidth',2)
    
    box on
    grid on
    
    %plot(K,z,'og','MarkerSize',10)

    %legend(k_names)
    
    %plot([10^-6,10^-2],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')

   
    xlabel(strcat('DPP \it K', '\rm (m/s)'))
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
   
    set(gca,'XScale','log')
%    set(gca,'FontSize',14)
%    set(gca,'YTickLabel',[]);
% 
%     xh = get(gca,'xlabel');
%     p = get(xh,'position'); % get the current position property
%     p(2) = 2*p(2) ;        % double the distance, 
%                            % negative values put the label below the axis
%     set(xh,'position',p) 
        set(gca,'FontSize',14)
    ylabel('Depth (m)')
    xlim([10^-6, 10^-2])
    



end
legend(sites)
end




