% Compare SOE with SumEch to K_SOE 
clear

baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];

exp = 2;

for k = 1:length(siteList)

    site = siteList{k}; 
    siteName = site
    
    [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);

    C_VC = 4200;
    
    % load data file
    [d, Dk{k}, T2ML, phi, z, SumEch{k},K_SOE{k}, log10K, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2_Ksoe(nmrName); 

    SOEest_fromK{k} = sqrt(K_SOE{k} ./ C_VC);
    
    % do reverse calc as check
    estSOEK{k} = C_VC .* (SOEest_fromK{k}).^2;
    
    
    %%%%%%%%% Change variable to Sum of Echoes for the bootstrap 
    lt = log10(SOEest_fromK{k}); 
    %lt = log10(SumEch{k});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nboot =  2000; % number of bootstrap samples

    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, log10K], Nboot, exp);    % n is fixed
   
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot)
    mediann(k) = median(n_boot)

    b_boot_all{k} = b_boot;
    
    %SOE_K_btstrp = medianb(k)*(SOEest_fromK{k}).^mediann(k);
    SOE_K_btstrp = medianb(k)*(SumEch{k}).^mediann(k);

    SumEch_k_estimates{k} = SOE_K_btstrp;
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate_btstrp(k) = median(estimateKdiffFactor(Dk{k}, SOE_K_btstrp, 1));
    errorEstimate_KSOE(k) = median(estimateKdiffFactor(Dk{k}, K_SOE{k}, 1));

end
