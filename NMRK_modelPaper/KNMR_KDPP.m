% Make NMR K vs  DPP K plot
clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants to use

%12/7/20 
%Using m = 1, n = 2 IF THIS IS CHANGED CHECK MAIN MODEL FILE CODE, NEED TO
%UPDATE VALUES IN SCRIPTS

m = 1;
n = 2;
siteList = [{"Site1-WellG5"},{"Site1-WellG6"},{"Site2-WellPN1"},{"Site2-WellPN2"}];

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
 
end
%%
% SOE Model

% ! SPECIFY PATH IN THIS FILE
computeSOE

SOEn = mediann;
SOEb = medianb;
SOEK = k_estimates;


%%
for k = 1:length(siteList)
   
   [SOE_errorSign,SOE_errorFactor] = estimateKdiffFactor_withSign(DPP_K{:,:,k},SOEK{:,k},1);
   SOE_errorSigns{:,k} = SOE_errorSign;
   SOE_errorFactors{:,k} = SOE_errorFactor;
   SOE_site{:,k} = repmat(siteList{k},length(SOE_errorSigns{:,k}),1);

end

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

Seeversm = [m m m m]';
Seeversn = [n n n n]';

for k = 1:length(siteList)
    
    site = siteList{k}
    [K,z,T2dist,T2logbins,SeeversK{k},k_mcmc,k_direct,bestFitMatrix,b_boot,totalErrorEstimate] = computeSeevers(site,Seeversn(k),Seeversm(k),figureson,wDirect);
   
   [Seevers_errorSign,Seevers_errorFactor] = estimateKdiffFactor_withSign(DPP_K{:,:,k},SeeversK{k},1);
   Seevers_errorSigns{:,k} = Seevers_errorSign;
   Seevers_errorFactors{:,k} = Seevers_errorFactor;
   Seevers_site{:,k} = repmat(siteList{k},length(SeeversK{k}),1);

end

%%
%TC model

TCm = [m m m m];
TCn = [n n n n];

% ! SPECIFY PATH IN THIS FILE
TC_gridSearch_wisc;

% this is set to 33 ms as of 7/8
cutoff = [cutoff cutoff cutoff cutoff];
cTC = totalcMatrix;

for k = 1:length(siteList)

    [TC_errorSign,TC_errorFactor] = estimateKdiffFactor_withSign(DPP_K{k},totalkTC{k}{:},1);
    TC_errorSigns{:,k} = TC_errorSign;
    TC_errorFactors{:,k} = TC_errorFactor;
    TC_site{:,k} = repmat(siteList{k},length(TC_errorSigns{:,k}),1);

end


%%
totalSDRK = vertcat(SDRK{:});
totalSeeversK = vertcat(SeeversK{:});
totalKGMK = vertcat(KGMK{:});
totalSOEK = vertcat(SOEK{:});
totalTCK = vertcat(totalkTC{:});
totalTCK = vertcat(totalTCK{:});

totalDPPK = vertcat(DPP_K{:});

totalDPPsites = vertcat(SDR_site{:});
totalSDRSites = vertcat(SDR_site{:});
[sortedDPPK, sortInd] = sort(totalDPPK);
sortedSDRSites_byDPP = totalSDRSites(sortInd);
sortedDPPSites_byDPP = totalDPPsites(sortInd);
sortedSDRK_byDPP = totalSDRK(sortInd);

figure()

hold on
    
grid on 
box on

scatter(totalDPPK, totalSDRK,80,'LineWidth',1.2)
scatter(totalDPPK, totalSeeversK,60,'+','LineWidth',1.2)
scatter(totalDPPK, totalKGMK,40,'^','LineWidth',1.2)
scatter(totalDPPK, totalSOEK,20,'<','LineWidth',1.2)
scatter(totalDPPK, totalTCK,30,'s','LineWidth',1.2)


plot(totalDPPK,totalDPPK,'k','LineWidth',2,'HandleVisibility','off')
plot(totalDPPK,totalDPPK*10,'k--','HandleVisibility','off')
plot(totalDPPK,totalDPPK*0.1,'k--','HandleVisibility','off')
plot([5*10^-5, 5*10^-5],[5*10^-8,1*10^0],'LineWidth',2)

legend('SDR','Seevers','KG','SOE','TC','Location','northwest')
ylim([4*10^-7,1*10^-2])
xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')
