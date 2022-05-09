% Make NMR-K grid plot 

close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants to use


%12/7/20 
%Using m = 1, n = 2 IF THIS IS CHANGED CHECK MAIN MODEL FILE CODE, NEED TO
%UPDATE VALUES IN SCRIPTS

m = 1;
n = 2;
siteList = [{"Site1-WellG5"},{"Site1-WellG6"},{"Site2-WellPN1"},{"Site2-WellPN2"}];
figureson = 0;

SDR_Kall = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_Kall = @(b,n,SOE) b.*(SOE).^n;
Seevers_Kall = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

SDRrange = [0.0132, 0.0213];
SOErange = [0.0077,0.016];
Seeversrange = [0.0106,0.0171];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% SDR model

% ! SPECIFY PATH IN THIS FILEs
runModelPairs

SDRm = [m m m m]';
SDRn = [n n n n]';
SDRb = squeeze(totalbMatrix(1,1,:));
SDRK = totalKMatrix;
DPP_K = DPP_KMatrix;

for k = 1:length(siteList)
  
   [SDR_errorSign,SDR_errorFactor] = estimateKdiffFactor_withSign(DPP_K{:,:,k},SDRK{:,:,k},1);
   SDR_site{:,k} = repmat(string(siteList{k}),length(SDR_errorSign),1);
   SDR_errorSigns{:,k} = SDR_errorSign;
   SDR_errorFactors{:,k} = SDR_errorFactor;
   
   SDRminK{k} = SDR_Kall(SDRrange(1),m,n,phi{k},T2ML{k});
   SDRmaxK{k} = SDR_Kall(SDRrange(2),m,n,phi{k},T2ML{k});
end

% SDR_T2ML = vertcat(T2ML{:});
% SDR_phi = vertcat(phi{:});
% 
% SDRminK = SDR_Kall(SDRrange(1),m,n,SDR_phi,SDR_T2ML);
% SDRmaxK = SDR_Kall(SDRrange(2),m,n,SDR_phi,SDR_T2ML);
%%
% SOE Model

% NOTE: FROM VISTA-CLARA, SHOULD USE SOE FROM K (KSOE) TO CALC SOE TO
% INCLUDE 500 MS WINDOWING. FOR THE WISCONSIN DATA, KSOE^2 = SUMECH, SO
% SUMECH IS ALREADY SQUARED!!
n = 1

% ! SPECIFY PATH IN THIS FILE
computeSOE

% n = 1 here is really n = 2!!!

SOEn = mediann;
SOEb = medianb;
SOEK = k_estimates;


%%
for k = 1:length(siteList)
   
    SOEminK{k} = SOE_Kall(SOErange(1),n,SumEch{k});
    SOEmaxK{k} = SOE_Kall(SOErange(2),n,SumEch{k});
 
    
   [SOE_errorSign,SOE_errorFactor] = estimateKdiffFactor_withSign(DPP_K{:,:,k},SOEK{:,k},1);
   SOE_errorSigns{:,k} = SOE_errorSign;
   SOE_errorFactors{:,k} = SOE_errorFactor;
   SOE_site{:,k} = repmat(siteList{k},length(SOE_errorSigns{:,k}),1);

end

n=2
%%
% KGM Model

% ! SPECIFY PATH IN THIS FILE
calcKGM

KGMK = KGM_K;
KGM_errorSigns = errorSign;
KGM_errorFactors = errorFactor;

for k =1:length(sites)
       KGM_sites{:,k} = repmat(sites{k},length(KGM_errorSigns{k}),1);
end
%%
% Seevers Model 

clear T2ML T2B phi

Seeversm = [m m m m]';
Seeversn = [n n n n]';
wDirect = 0;

for k = 1:length(siteList)
    
    site = char(siteList{k})
    [K,z,T2dist,T2logbins,SeeversK{k},k_mcmc,k_direct,bestFitMatrix,b_boot,idealbs,totalErrorEstimate,T2ML{k},T2B{k},phi{k}] = computeSeevers(site,Seeversn(k),Seeversm(k),figureson,wDirect);
   
   [Seevers_errorSign,Seevers_errorFactor] = estimateKdiffFactor_withSign(DPP_K{k},SeeversK{k},1);
   Seevers_errorSigns{:,k} = Seevers_errorSign;
   Seevers_errorFactors{:,k} = Seevers_errorFactor;
   Seevers_site{:,k} = repmat(siteList{k},length(SeeversK{k}),1);
   
   SeeversminK{k} = Seevers_Kall(Seeversrange(1),m,n,T2ML{k},T2B{k},phi{k});
   SeeversmaxK{k} = Seevers_Kall(Seeversrange(2),m,n,T2ML{k},T2B{k},phi{k});

end



%%

figure(1)
t = tiledlayout(4,4);
columnNames = [{"SDR"},{"Seevers"},{"SOE"},{"KG"}];

minDepth = 1;
maxDepth = 18;

DPP_depths = finalDepths;
DPP_depthsAll = vertcat(DPP_depths{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1
    

nexttile
hold on
grid on 
box on
grid minor

scatter(SDRminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SDRmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SDRK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
       
ylabel('Depth (m)') 

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
title(columnNames{1})

set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

nexttile
hold on
grid on 
box on
grid minor

scatter(SeeversminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SeeversmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SeeversK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)

set(gca,'YTickLabel',[]);

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
title(columnNames{2})
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);


nexttile
hold on
grid on 
box on
grid minor

scatter(SOEminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SOEmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SOEK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)

set(gca,'YTickLabel',[]);

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
title(columnNames{3})
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

nexttile
hold on
grid on 
box on
grid minor

scatter(KGMK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])
scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
set(gca,'YTickLabel',[]);


%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
title(columnNames{4})
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2
    

nexttile
hold on
grid on 
box on
grid minor

scatter(SDRminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SDRmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SDRK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
       
ylabel('Depth (m)') 

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])    

set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

nexttile
hold on
grid on 
box on
grid minor

scatter(SeeversminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SeeversmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SeeversK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)

set(gca,'YTickLabel',[]);

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

nexttile
hold on
grid on 
box on
grid minor

scatter(SOEminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SOEmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SOEK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)

set(gca,'YTickLabel',[]);

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

nexttile
hold on
grid on 
box on
grid minor

scatter(KGMK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])
scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
set(gca,'YTickLabel',[]);


%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 3
    

nexttile
hold on
grid on 
box on
grid minor

scatter(SDRminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SDRmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SDRK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
       
ylabel('Depth (m)') 

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])    

set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

nexttile
hold on
grid on 
box on
grid minor

scatter(SeeversminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SeeversmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SeeversK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)

set(gca,'YTickLabel',[]);

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

nexttile
hold on
grid on 
box on
grid minor

scatter(SOEminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SOEmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SOEK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)

set(gca,'YTickLabel',[]);

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

nexttile
hold on
grid on 
box on
grid minor

scatter(KGMK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])
scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
set(gca,'YTickLabel',[]);


%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])
set(gca,'XTickLabel',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 4
    

nexttile
hold on
grid on 
box on
grid minor

scatter(SDRminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SDRmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SDRK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
       
ylabel('Depth (m)') 

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])    

set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])


nexttile
hold on
grid on 
box on
grid minor

scatter(SeeversminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SeeversmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SeeversK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)

set(gca,'YTickLabel',[]);

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])

nexttile
hold on
grid on 
box on
grid minor

scatter(SOEminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
scatter(SOEmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)

scatter(SOEK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])

scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)

set(gca,'YTickLabel',[]);

%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])

nexttile
hold on
grid on 
box on
grid minor

scatter(KGMK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])
scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
set(gca,'YTickLabel',[]);


%xlabel(strcat('DPP \it K', '\rm (m/s)'))
set(gca, 'YDir','reverse')
ylim([minDepth,maxDepth])
set(gca,'XScale','log')
xlim([3*10^-6, 2*10^-3])

% for k = 1:length(siteList)
%     
%     nexttile
%     hold on
%     grid on 
%     box on
%     grid minor
% 
%     
%     scatter(SeeversminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
%     scatter(SeeversmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
% 
%     scatter(SeeversK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])
% 
%     scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
%     
% 
%    
%     set(gca, 'YDir','reverse')
%     ylim([minDepth,maxDepth])
%    
%     set(gca,'XScale','log')
% 
%     xlim([3*10^-6, 2*10^-3])
%     
%     if k == 1
%        ylabel('Depth (m)') 
%     else
%        set(gca,'YTickLabel',[]);
%     end
% 
%      set(gca,'XTickLabel',[]);
%     
% end
% 
% for k = 1:length(siteList)
%     
%     nexttile
%     hold on
%     grid on 
%     box on
%     grid minor
% 
%     scatter(SOEminK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
%     scatter(SOEmaxK{k},DPP_depths{k},90,[0.5 0.5 0.5],'|','LineWidth',1.2)
% 
%     scatter(SOEK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])
% 
%     scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
%    
%     %xlabel(strcat('DPP \it K', '\rm (m/s)'))
%     set(gca, 'YDir','reverse')
%     ylim([minDepth,maxDepth])
%    
%     set(gca,'XScale','log')
% 
%     xlim([3*10^-6, 2*10^-3])
%     
%     if k == 1
%        ylabel('Depth (m)') 
%     else
%        set(gca,'YTickLabel',[]);
%     end
%     set(gca,'XTickLabel',[]);
% 
%     
% end
% 
% for k = 1:length(siteList)
%     
%     nexttile
%     hold on
%     grid on 
%     grid minor
%     box on
%     
%     scatter(KGMK{k},DPP_depths{k},20,'Filled','MarkerEdgeColor', [0, 0.4470, 0.7410],'MarkerFaceColor', [0, 0.4470, 0.7410])
%     scatter(DPP_K{k}, DPP_depths{k},30,'+','r','LineWidth',1.2)
%    
%     %xlabel(strcat('DPP \it K', '\rm (m/s)'))
%     set(gca, 'YDir','reverse')
%     ylim([minDepth,maxDepth])
%    
%     set(gca,'XScale','log')
% 
%     xlim([3*10^-6, 2*10^-3])
%     
%     if k == 1
%        ylabel('Depth (m)') 
%     else
%        set(gca,'YTickLabel',[]);
%     end
%     
%     xlabel('K (m/s)')
% 
%     
% end

