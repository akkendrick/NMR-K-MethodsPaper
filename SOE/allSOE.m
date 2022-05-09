baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

alln = [1 2];
siteList = [{'Site1-WellG5'} {'Site1-WellG6'}  {'Site2-WellPN1'} {'Site2-WellPN2'}];
for j = 1:length(alln)
    disp('Current n')
    currentn = alln(j)
    
    for k = 1:length(siteList)

        site = siteList{k}; 
        siteName = site

        [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);

        % load data file
        [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
            SumEch_twm_3s] = loadnmrdata2(nmrName); 

        logSumEch = log10(SumEch); 

        %%%%%%%%% Change T2 variable to Sum of Echoes for the inversions. 
        lt = logSumEch; 
        T2ML = SumEch; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Nboot =  2000; % number of bootstrap samples

        %[b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
        [b_boot, n_boot] = bootstrap_fun([lt, kk], Nboot, currentn);    % this is
        % configured such that n is fixed!
        

        % Need to fix n here too!
        meanb{j,k} = mean(b_boot); 
        meann{j,k} = currentn; 

        medianb{j,k} = median(b_boot)
        mediann{j,k} = currentn;

        b_boot_all{j,k} = b_boot;

        SOE_K{j,k} = medianb{j,k}*(SumEch).^mediann{j,k};
        k_estimates{j,k} = SOE_K{j,k};
        
        DPP_K{j,k} = Dk;

        %plotKwithDepth(Dk,z,T2dist,T2logbins,SOE_K,{'DPP','SOE'},{'+'})
        errorEstimate{j,k} = sum(computeError(Dk, SOE_K{j,k}));

    end

    allSOEK{j,:} = vertcat(k_estimates{j,:})
    allSOEError{j,:} = vertcat(errorEstimate{j,:})
    allDPPK{j,:} = vertcat(DPP_K{j,:})
end

%% Plot data
figure(1)

hold on
    
grid on
box on

scatter(allDPPK{1,:}, allSOEK{1,:},60,'Filled')
scatter(allDPPK{1,:}, allSOEK{2,:},60,'Filled')

legend('SOE n = 1','SOE n = 2','Location','northwest')

plot(allDPPK{1,:},allDPPK{1,:},'k','LineWidth',2,'HandleVisibility','off')
plot(allDPPK{1,:},allDPPK{1,:}*10,'k:','HandleVisibility','off')
plot(allDPPK{1,:},allDPPK{1,:}*0.1,'k:','HandleVisibility','off')

ylim([5*10^-8,1*10^0])
xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')

%% Compute Avg K diff factor
disp('SOE n = 1 mean K diff factor')
mean(vertcat(errorEstimate{1,:}))

disp('SOE n = 2 mean K diff factor')
mean(vertcat(errorEstimate{2,:}))
