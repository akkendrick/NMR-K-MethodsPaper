% Plot all different types of SOE

% Compare SOE with SumEch to K_SOE 
clear

baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];

n = 1;

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
    %lt = log10(SOEest_fromK{k}); 
    lt = log10(SumEch{k});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nboot =  2000; % number of bootstrap samples

    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, log10K], Nboot, n);    % n is fixed
   
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot)
    mediann(k) = median(n_boot)

    b_boot_all{k} = b_boot;
    
    %SOE_K_btstrp = medianb(k)*(SOEest_fromK{k}).^mediann(k);
    SOE_K_btstrp_n1 = medianb(k)*(SumEch{k}).^mediann(k);

    SumEch_k_estimates_n1{k} = SOE_K_btstrp_n1;
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate_btstrp_n1(k) = median(estimateKdiffFactor(Dk{k}, SOE_K_btstrp_n1, 1));
    errorEstimate_KSOE(k) = median(estimateKdiffFactor(Dk{k}, K_SOE{k}, 1));

end

n = 1;

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nboot =  2000; % number of bootstrap samples

    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, log10K], Nboot, n);    % n is fixed
   
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot)
    mediann(k) = median(n_boot)

    b_boot_all{k} = b_boot;
    
    %SOE_K_btstrp = medianb(k)*(SOEest_fromK{k}).^mediann(k);
    SOEest_k_est = medianb(k)*(SumEch{k}).^mediann(k);

    SOEest_k_estimates_n1{k} = SOEest_k_est;
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate_SOEest_btstrp_n1(k) = median(estimateKdiffFactor(Dk{k}, SOEest_k_est, 1));
    errorEstimate_KSOE(k) = median(estimateKdiffFactor(Dk{k}, K_SOE{k}, 1));

end
%%
n = 2;

for k = 1:length(siteList)

    site = siteList{k}; 
    siteName = site
    
    [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);

    C_VC = 4200;
    
    % load data file
    [d, Dk{k}, T2ML, phi, z, SumEch{k},K_SOE{k}, log10K, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2_Ksoe(nmrName); 

    DPP_K{k} = Dk{k}
    
    SOEest_fromK{k} = sqrt(K_SOE{k} ./ C_VC);
    
    % do reverse calc as check
    estSOEK{k} = C_VC .* (SOEest_fromK{k}).^2;
    
    
    %%%%%%%%% Change variable to Sum of Echoes for the bootstrap 
    %lt = log10(SOEest_fromK{k}); 
    lt = log10(SumEch{k});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nboot =  2000; % number of bootstrap samples

    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, log10K], Nboot, n);    % n is fixed
   
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot)
    mediann(k) = median(n_boot)

    b_boot_all{k} = b_boot;
    
    %SOE_K_btstrp = medianb(k)*(SOEest_fromK{k}).^mediann(k);
    SOE_K_btstrp = medianb(k)*(SumEch{k}).^mediann(k);

    SumEch_k_estimates_n2{k} = SOE_K_btstrp;
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate_btstrp_n2(k) = median(estimateKdiffFactor(Dk{k}, SOE_K_btstrp, 1));
    errorEstimate_KSOE(k) = median(estimateKdiffFactor(Dk{k}, K_SOE{k}, 1));

end

n=2;

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nboot =  2000; % number of bootstrap samples

    %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    [b_boot, n_boot] = bootstrap_fun([lt, log10K], Nboot, n);    % n is fixed
   
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 
    
    medianb(k) = median(b_boot)
    mediann(k) = median(n_boot)

    b_boot_all{k} = b_boot;
    
    %SOE_K_btstrp = medianb(k)*(SOEest_fromK{k}).^mediann(k);
    SOEest_k_est = medianb(k)*(SumEch{k}).^mediann(k);

    SOEest_k_estimates_n2{k} = SOEest_k_est;
    
    %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
    errorEstimate_SOEest_btstrp_n2(k) = median(estimateKdiffFactor(Dk{k}, SOEest_k_est, 1));
    errorEstimate_KSOE(k) = median(estimateKdiffFactor(Dk{k}, K_SOE{k}, 1));

end
%%
totalDPPK = vertcat(DPP_K{:});
totalSOEK_VC_n2 = vertcat(K_SOE{:});
totalSOEK_SumEch_n2 = vertcat(SumEch_k_estimates_n2{:});
totalSOEK_SumEch_n1 = vertcat(SumEch_k_estimates_n1{:});
totalSOEK_SOEest_n2 = vertcat(SOEest_k_estimates_n2{:});
totalSOEK_SOEest_n1 = vertcat(SOEest_k_estimates_n1{:});

figure()

hold on
    
grid on 
box on

scatter(totalDPPK, totalSOEK_VC_n2,60,'Filled')
scatter(totalDPPK, totalSOEK_SumEch_n2,60,'Filled')
scatter(totalDPPK, totalSOEK_SumEch_n1,60,'Filled')
scatter(totalDPPK, totalSOEK_SOEest_n2,60,'Filled')
scatter(totalDPPK, totalSOEK_SOEest_n1,60,'Filled')

plot(totalDPPK,totalDPPK,'k','LineWidth',2,'HandleVisibility','off')
plot(totalDPPK,totalDPPK*10,'k:','HandleVisibility','off')
plot(totalDPPK,totalDPPK*0.1,'k:','HandleVisibility','off')
plot([5*10^-5, 5*10^-5],[5*10^-8,1*10^0],'LineWidth',2)

legend('Vista Clara KSOE, n=2, Window=500ms','SumEch SOE, n=2, Window=???','SumEch SOE, n=1, Window=???',...
    'K SOE Est, n=2, Window=500ms','K SOE Est, n=1, Window=500ms','Location','northwest')
ylim([5*10^-8,1*10^0])
xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')