% computeSOE
clear

% Analyze SOE by plotting K_SOE vs K_DPP
%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

n = 2;

siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];

for i = 1:length(siteList)
    siteName = siteList{i};
    [T2dist{i}, T2logbins{i}, nmrName{i}] = loadRawNMRdata(siteName);
    
    [d{i}, K{i}, T2ML{i}, phi{i}, z{i}, SumEch{i},K_SOE{i}, logK{i}, logT2ML{i}, logPhi{i}, SumEch_3s{i}, SumEch_twm{i}, ...
    SumEch_twm_3s{i}] = loadnmrdata2_Ksoe(nmrName{i}); 
    
    depthsAll = z{i};
    finalDepths = depthsAll;

if (siteName == "Site1-WellG6")


        depthCutoff = 5.8;

        K{i} = K{i}(depthsAll>depthCutoff);
        phi{i} = phi{i}(depthsAll>depthCutoff);
        T2ML{i} = T2ML{i}(depthsAll>depthCutoff);
        SumEch{i} = SumEch{i}(depthsAll>depthCutoff);
        logK{i} = logK{i}(depthsAll>depthCutoff);
        logT2ML{i} = logT2ML{i}(depthsAll>depthCutoff);
        logPhi{i} = logPhi{i}(depthsAll>depthCutoff);
        T2dist{i} = T2dist{i}(T2dist{i}(:,1)>depthCutoff,:);
        K_SOE{i} = K_SOE{i}(depthsAll>depthCutoff);

    elseif (siteName == "Site1-WellG5")


        depthCutoff = 4;

        K{i} = K{i}(depthsAll>depthCutoff);
        phi{i} = phi{i}(depthsAll>depthCutoff);
        T2ML{i} = T2ML{i}(depthsAll>depthCutoff);
        SumEch{i} = SumEch{i}(depthsAll>depthCutoff);
        logK{i} = logK{i}(depthsAll>depthCutoff);
        logT2ML{i} = logT2ML{i}(depthsAll>depthCutoff);
        logPhi{i} = logPhi{i}(depthsAll>depthCutoff);
        T2dist{i} = T2dist{i}(T2dist{i}(:,1)>depthCutoff,:);
        K_SOE{i} = K_SOE{i}(depthsAll>depthCutoff);

    end
end

Kall = vertcat(K{:});
logK_all = vertcat(logK{:});
logT2ML_all = vertcat(logT2ML{:});
logPhi_all = vertcat(logPhi{:});
T2logbins_all = T2logbins{1};
T2dist_all = vertcat(T2dist{:});
SumEch_all = vertcat(SumEch{:});
Zall = vertcat(z{:});
K_SOE_all = vertcat(K_SOE{:})

logSumEch = log10(SumEch_all); 
logK = log10(Kall);

C_VC = 4200;
SOEest_fromK = sqrt(K_SOE_all ./ C_VC);


%%%%%%%%% Change variable to Sum of Echoes for the inversions. 
lt = log10(SOEest_fromK); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nboot =  2000; % number of bootstrap samples

[b_boot, n_boot] = bootstrap_fun([lt, logK], Nboot, n);    % n is fixed

disp('Median prefactor for SOE is:')
medianb = median(b_boot)

SOE_K = medianb*(SOEest_fromK).^n;
k_estimates = SOE_K;

disp('Median K diff factor for SOE is:')
errorEstimate = median(estimateKdiffFactor(Kall, SOE_K,1))
    
