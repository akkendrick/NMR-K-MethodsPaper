% computeSOE
%clear

% Analyze SOE by plotting K_SOE vs K_DPP
%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
baseDir = 'C:\Users\kenta\Dropbox\Research\Alex-Rosemary\Papers\Kmodel_comparision_paper\howTo\Common\Field Data\USGS Data\';

n = 1;

k_estimates = [];
k_names = {'1:1','SOE n=2','SOE n=1','SOE n=0'};
k_sym = {'+','*','o'};

siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];

% siteList = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%   'dpnmr_leque_east','dpnmr_leque_west'};


for k = 1:length(siteList)

    site = siteList{k}; 
    siteName = site
    
    [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);

    % load data file
    [d, DK{k}, T2ML, phi, z, SumEch{k},K_SOE{k}, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2_Ksoe(nmrName); 
    
    depthsAll = z;
    
% Filter data by depth
if (site == "Site1-WellG6")


    depthCutoff = 5.8;

    DK{k} = DK{k}(depthsAll>depthCutoff);
    phi = phi(depthsAll>depthCutoff);
    T2ML = T2ML(depthsAll>depthCutoff);
    SumEch{k} = SumEch{k}(depthsAll>depthCutoff);
    kk = kk(depthsAll>depthCutoff);

    
elseif (site == "Site1-WellG5")


    depthCutoff = 4;
    
    DK{k} = DK{k}(depthsAll>depthCutoff);
    phi = phi(depthsAll>depthCutoff);
    T2ML = T2ML(depthsAll>depthCutoff);
    SumEch{k} = SumEch{k}(depthsAll>depthCutoff);
    kk = kk(depthsAll>depthCutoff);


end

    logSumEch = log10(SumEch{k}); 

    %%%%%%%%% Change T2 variable to Sum of Echoes for the bootstrap. 
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
    
    SOE_K = medianb(k)*(SumEch{k}).^mediann(k);
    k_estimates{k} = SOE_K;
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate(k) = median(estimateKdiffFactor(DK{k}, SOE_K, 1));
    
end

save('SOE_dat_n1.mat','siteList','b_boot_all','medianb','mediann','k_estimates','errorEstimate')
% plotKestKdpp(Dk,k_estimates,k_names,k_sym)
% title(siteName)
