% Make prefactor plot

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants and equations to use

%12/7/20 
%Using m = 1, n = 2 IF THIS IS CHANGED CHECK MAIN MODEL FILE CODE, NEED TO
%UPDATE VALUES IN SCRIPTS

m = 1;
n = 2;
siteList = [{'Site1-WellG5'},{'Site1-WellG6'},{'Site2-WellPN1'},{'Site2-WellPN2'}];

% SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
% SOE_K = @(b,n,SOE) b.*(SOE).^n;
% Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% SDR model

% ! SPECIFY PATH IN THIS FILEs
runModelPairs

SDR_b = @(m,n,phi,T2ML,K) K ./ ((phi.^m).*(T2ML).^n);

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
   SDR_bs{:,k} = SDR_b(m,n,phi{k},T2ML{k},DPP_K{:,:,k});
 
end
%%
% SOE Model

% ! SPECIFY PATH IN THIS FILE
computeSOE

SOE_b = @(n,K,SOE) K./(SOE).^n;

SOEn = mediann;
SOEb = medianb;
SOEK = k_estimates;

for k = 1:length(siteList)
   
   [SOE_errorSign,SOE_errorFactor] = estimateKdiffFactor_withSign(DK{k},SOEK{k},1);
   SOE_errorSigns{:,k} = SOE_errorSign;
   SOE_errorFactors{:,k} = SOE_errorFactor;
   SOE_site{:,k} = repmat(siteList{k},length(SOE_errorSigns{:,k}),1);

   SOE_bs{:,k} = SOE_b(n,DK{k},SumEch{k});

end

%%
% Seevers Model 
m = 1;
n = 2;

Seeversm = [m m m m]';
Seeversn = [n n n n]';


% Specify options
figureson = 0;
wDirect = 0;

% Computing bs is in the computeSeevers file
for k = 1:length(siteList)
    
    site = siteList{k}
    
    [DPP_K{k},z,T2dist,T2logbins,SeeversK{k},k_mcmc,k_direct,...
        bestFitMatrix,b_boot,Seevers_bs{k},totalErrorEstimate] = computeSeevers(site,...
        Seeversn(k),Seeversm(k),figureson,wDirect);
   
   [Seevers_errorSign,Seevers_errorFactor] = estimateKdiffFactor_withSign(DPP_K{k},...
       SeeversK{k},1);
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
%

%%
totalDPPK = vertcat(DPP_K{:});


% Save all the bs
allTCb = vertcat(TC_bs{:});
allSDRb = vertcat(SDR_bs{:});
allSeeversb = vertcat(Seevers_bs{:});
allSOEb = vertcat(SOE_bs{:});

% Plot all the prefactors

figure()
 
hold on
    
grid on 
box on

scatter(totalDPPK, allSDRb,60,'Filled')
scatter(totalDPPK, allSeeversb,60,'Filled')
scatter(totalDPPK, allSOEb,60,'Filled')
scatter(totalDPPK, allTCb,60,'Filled')

legend('SDR','Seevers','SOE','TC 33 ms cutoff','Location','northwest')
ylim([5*10^-8,1*10^0])
xlabel('DPP K (m/s)')
ylabel('Calibration Prefactors') 

set(gca,'FontSize',16)
% set(gca,'XScale','log')
% set(gca,'YScale','log')

%%