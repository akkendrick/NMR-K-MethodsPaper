% Run Timur-Coates Model Estimates
% Range over pairs of m and n values
 close all
 clear

load enso 

siteList = [{'Site1-WellG5'} {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

TC_b = @(K,phi,frac) K ./ ((phi).*(frac).^2);

cutoff = 33*10^-3;

m = 1;
n = 2;

saveData = 0;

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
T2logbins_all = T2logbins{1};
T2dist_all = vertcat(T2dist{:});
Zall = vertcat(z{:});

[b_boot,k_boot,totalErrorEstimate] = computeTCperm_all(Kall,Zall,T2dist_all,...
     T2logbins_all,logPhi_all,n,m,cutoff);


    
    

