% Estimate KGM parameters from all wells
clear

baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

n = 2;
m = 1;
Tb= @(T) 3.3 + .044.*(T - 35);     

siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];

for i = 1:length(siteList)
    siteName = siteList{i};
    [T2dist{i}, T2logbins, nmrName{i}] = loadRawNMRdata(siteName);
    
    [d, DPP_K{i}, T2ML{i}, phi{i}, z{i}, SumEch{i},K_SOE{i}, logK{i}, logT2ML{i}, logPhi{i}, SumEch_3s{i}, SumEch_twm{i}, ...
    SumEch_twm_3s{i}] = loadnmrdata2_Ksoe(nmrName{i}); 

%     % Compute the maximum T2
% 
%     logT2peakVal = []
%     nmrDepths = d(:,1);
%     for j = 1:length(nmrDepths)
%         depth = nmrDepths(j)
%         T2depths = T2dist{i}(:,1);
%         T2ind = (round(T2depths,4) == depth)
%         T2interval = T2dist{i}(T2ind,:);
%         [T2peakAmp, T2peakInd] = max(T2interval(2:end));
%         logT2peakVal(j) = T2logbins(T2peakInd);
%     end
%     
%     logT2peakVal_step{i} = logT2peakVal';
%     T2peakVal{i} = 10.^(logT2peakVal');
    logK{i} = log10(DPP_K{i});

    depthsAll = z{i};

    
    % Filter data by depth
    if (siteName == "Site1-WellG6")


        depthCutoff = 5.8;

        DPP_K{i} = DPP_K{i}(depthsAll>depthCutoff);
        phi{i} = phi{i}(depthsAll>depthCutoff);
        T2ML{i} = T2ML{i}(depthsAll>depthCutoff);
        SumEch{i} = SumEch{i}(depthsAll>depthCutoff);
        logK{i} = logK{i}(depthsAll>depthCutoff);
        logT2ML{i} = logT2ML{i}(depthsAll>depthCutoff);
        logPhi{i} = logPhi{i}(depthsAll>depthCutoff);
        T2dist{i} = T2dist{i}(T2dist{i}(:,1)>depthCutoff,:);

    elseif (siteName == "Site1-WellG5")


        depthCutoff = 4;

        DPP_K{i} = DPP_K{i}(depthsAll>depthCutoff);
        phi{i} = phi{i}(depthsAll>depthCutoff);
        T2ML{i} = T2ML{i}(depthsAll>depthCutoff);
        SumEch{i} = SumEch{i}(depthsAll>depthCutoff);
        logK{i} = logK{i}(depthsAll>depthCutoff);
        logT2ML{i} = logT2ML{i}(depthsAll>depthCutoff);
        logPhi{i} = logPhi{i}(depthsAll>depthCutoff);
        T2dist{i} = T2dist{i}(T2dist{i}(:,1)>depthCutoff,:);
    end
    
    logK{i} = log10(DPP_K{i});
    %T2B_site{i} =  repmat(2.0680,length(DPP_K{i}),1);
    T2B_site{i} =  repmat(2.2440,length(DPP_K{i}),1);

end

T2B_all = vertcat(T2B_site{:});
logT2ML_all = log10(vertcat(T2ML{:}));
% logT2peakVal_all = vertcat(logT2peakVal_step{:});
phi_all = vertcat(phi{:});
logK_all = vertcat(logK{:});
DPP_K_all = vertcat(DPP_K{:});
    
[KGM_lk, bestTau, bestRho, r] = grid_search_kgm(logT2ML_all, phi_all,logK_all, m, T2B_all);

disp(strcat('Best tau for KGM is:', string(bestTau)))

disp(strcat('Best rho for KGM is:', string(bestRho)))

KGM_K = 10.^KGM_lk;

[errorSign, errorFactor] = estimateKdiffFactor_withSign(DPP_K_all,KGM_K,1);

disp(strcat('Median K diff factor for KGM is:', string(median(errorFactor))))
medianErrorFactor = median(errorFactor);


%% Plot results
totalKGMK = KGM_K;
totalDPPK = DPP_K_all;

figure(1)

hold on
    
grid on 
box on

scatter(totalDPPK, totalKGMK,60,'Filled')

plot(totalDPPK,totalDPPK,'k','LineWidth',2,'HandleVisibility','off')
plot(totalDPPK,totalDPPK*10,'k--','HandleVisibility','off')
plot(totalDPPK,totalDPPK*0.1,'k--','HandleVisibility','off')

legend('KGM','Location','northwest')
ylim([5*10^-8,1*10^0])
xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')