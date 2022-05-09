% Compare SOE with SumEch to K_SOE 
clear

baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];

n = 1;

for k = 1:length(siteList)

    site = siteList{k}; 
    siteName = site
    
    [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);

    % load data file
    [d, Dk{k}, T2ML, phi, z, SumEch{k},K_SOE{k}, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2_Ksoe(nmrName); 

    logSumEch = log10(SumEch{k}); 

    %%%%%%%%% Change T2 variable to Sum of Echoes for the inversions. 
    lt = logSumEch; 
    T2ML = SumEch{k}; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nboot =  2000; % number of bootstrap samples

    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, kk], Nboot, n);    % n is fixed
   
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot)
    mediann(k) = median(n_boot)

    b_boot_all{k} = b_boot;
    
    SumEch_SOE_K = medianb(k)*(SumEch{k}).^mediann(k);
    SumEch_k_estimates{k} = SumEch_SOE_K;
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate_SumEch(k) = sum(estimateKdiffFactor(Dk{k}, SumEch_SOE_K, 1));
    errorEstimate_KSOE(k) = sum(estimateKdiffFactor(Dk{k}, K_SOE{k}, 1));

end
