% Run Seevers Model for G6


close all
clear

siteList = [{'Site1-WellG6'}];

m = [1];
n = [2];
depthCutoff = 5.8; % depth that separates layers in m

for i = 1:length(siteList)
    siteName = siteList{i};
    [T2dist{i}, T2logbins{i}, nmrName{i}] = loadRawNMRdata(siteName);
    [d{i}, K{i}, T2ML{i}, phi{i}, z{i}, SumEch{i}, logK{i}, logT2ML{i}, logPhi{i}, SumEch_3s{i}, SumEch_twm{i}, ...
    SumEch_twm_3s{i}] = loadnmrdata2(nmrName{i}); 
end

Kall = vertcat(K{:});
logK_all = vertcat(logK{:});
logT2ML_all = vertcat(logT2ML{:});
logPhi_all = vertcat(logPhi{:});
depthsAll = vertcat(z{i});

[b_boot,k_boot,totalErrorEstimate] = computeSeevers_G6(Kall,logK_all,...
    logT2ML_all,logPhi_all,n,m,depthsAll,depthCutoff);

%%
figure()

hold on
    
grid on 
box on

scatter(Kall, k_boot,60,'Filled')

plot(Kall,Kall,'k','LineWidth',2,'HandleVisibility','off')
plot(Kall,Kall*10,'k--','HandleVisibility','off')
plot(Kall,Kall*0.1,'k--','HandleVisibility','off')

legend('Seevers','Location','northwest')
ylim([5*10^-8,1*10^0])
xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')