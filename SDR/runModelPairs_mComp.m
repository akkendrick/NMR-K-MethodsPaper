% Compare SDR m = 0 and m = 1
clear
close all

%% m = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wisconsin Data
siteList = [{'Site1-WellG5'}...
   {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fitting variables and K range
m = 0;
n = 2;
Kcutoff = 10^-9;
Nboot = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
for i = 1:length(siteList)
    siteName = siteList{i};
    disp(strcat('On site: ', siteName))
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(siteName);

    [d{i}, K{i}, T2ML{i}, phi{i}, z{i}, SumEch{i}, logK{i}, logT2ML{i}, logPhi{i}, SumEch_3s{i}, SumEch_twm{i}, ...
    SumEch_twm_3s{i}] = loadnmrdata2(nmrName); 

    T2depths_site{i} = T2dist(:,1);
    T2data_site{i} = T2dist(:,2:end);
    T2dist_site{i} = T2dist;

end

allK = vertcat(K{:});
allT2ML = vertcat(T2ML{:});
allPhi = vertcat(phi{:});

alllogK = vertcat(logK{:});
alllogT2ML = vertcat(logT2ML{:});
alllogPhi = vertcat(logPhi{:});

% Change this line to do > or <
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kind = allK > Kcutoff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allK_filt = allK(Kind);
allT2ML_filt = allT2ML(Kind);
allPhi_filt = allPhi(Kind);

alllogK_filt = alllogK(Kind);
alllogT2ML_filt = alllogT2ML(Kind);
alllogPhi_filt = alllogPhi(Kind);

    
[b_boot, n_boot, m_boot] = bootstrap_fixed_m_n([alllogT2ML_filt, alllogPhi_filt, alllogK_filt], Nboot, n, m);   % m, n fixed

median_b_m0 = median(b_boot)
NMR_K_m0 = median_b_m0*(allPhi_filt.^m).*(allT2ML_filt).^n;

[NMR_K_m0_errorSign,NMR_K_m0_errorFactor] = estimateKdiffFactor_withSign(allK,NMR_K_m0,1);

%% m = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wisconsin Data
siteList = [{'Site1-WellG5'}...
   {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fitting variables and K range
m = 1;
n = 2;
Kcutoff = 10^-9;
Nboot = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
for i = 1:length(siteList)
    siteName = siteList{i};
    disp(strcat('On site: ', siteName))
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(siteName);

    [d{i}, K{i}, T2ML{i}, phi{i}, z{i}, SumEch{i}, logK{i}, logT2ML{i}, logPhi{i}, SumEch_3s{i}, SumEch_twm{i}, ...
    SumEch_twm_3s{i}] = loadnmrdata2(nmrName); 

    T2depths_site{i} = T2dist(:,1);
    T2data_site{i} = T2dist(:,2:end);
    T2dist_site{i} = T2dist;

end

allK = vertcat(K{:});
allT2ML = vertcat(T2ML{:});
allPhi = vertcat(phi{:});

alllogK = vertcat(logK{:});
alllogT2ML = vertcat(logT2ML{:});
alllogPhi = vertcat(logPhi{:});

% Change this line to do > or <
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kind = allK > Kcutoff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allK_filt = allK(Kind);
allT2ML_filt = allT2ML(Kind);
allPhi_filt = allPhi(Kind);

alllogK_filt = alllogK(Kind);
alllogT2ML_filt = alllogT2ML(Kind);
alllogPhi_filt = alllogPhi(Kind);

    
[b_boot, n_boot, m_boot] = bootstrap_fixed_m_n([alllogT2ML_filt, alllogPhi_filt, alllogK_filt], Nboot, n, m);   % m, n fixed

median_b_m1 = median(b_boot)
NMR_K_m1 = median_b_m1*(allPhi_filt.^m).*(allT2ML_filt).^n;
[NMR_K_m1_errorSign,NMR_K_m1_errorFactor] = estimateKdiffFactor_withSign(allK,NMR_K_m1,1);


%% m = 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wisconsin Data
siteList = [{'Site1-WellG5'}...
   {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fitting variables and K range
m = 4;
n = 2;
Kcutoff = 10^-9;
Nboot = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
for i = 1:length(siteList)
    siteName = siteList{i};
    disp(strcat('On site: ', siteName))
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(siteName);

    [d{i}, K{i}, T2ML{i}, phi{i}, z{i}, SumEch{i}, logK{i}, logT2ML{i}, logPhi{i}, SumEch_3s{i}, SumEch_twm{i}, ...
    SumEch_twm_3s{i}] = loadnmrdata2(nmrName); 

    T2depths_site{i} = T2dist(:,1);
    T2data_site{i} = T2dist(:,2:end);
    T2dist_site{i} = T2dist;

end

allK = vertcat(K{:});
allT2ML = vertcat(T2ML{:});
allPhi = vertcat(phi{:});

alllogK = vertcat(logK{:});
alllogT2ML = vertcat(logT2ML{:});
alllogPhi = vertcat(logPhi{:});

% Change this line to do > or <
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kind = allK > Kcutoff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allK_filt = allK(Kind);
allT2ML_filt = allT2ML(Kind);
allPhi_filt = allPhi(Kind);

alllogK_filt = alllogK(Kind);
alllogT2ML_filt = alllogT2ML(Kind);
alllogPhi_filt = alllogPhi(Kind);

    
[b_boot, n_boot, m_boot] = bootstrap_fixed_m_n([alllogT2ML_filt, alllogPhi_filt, alllogK_filt], Nboot, n, m);   % m, n fixed

median_b_m4 = median(b_boot)
NMR_K_m4 = median_b_m4*(allPhi_filt.^m).*(allT2ML_filt).^n;

[NMR_K_m4_errorSign,NMR_K_m4_errorFactor] = estimateKdiffFactor_withSign(allK,NMR_K_m4,1);
%% Now plot results
totalDPPK = allK;

 figure()

hold on
    
grid on 
box on


scatter(totalDPPK, NMR_K_m0,60,'Filled')
scatter(totalDPPK, NMR_K_m1,60,'Filled')
scatter(totalDPPK, NMR_K_m4,60,'Filled')

plot(totalDPPK,totalDPPK,'k','LineWidth',2,'HandleVisibility','off')
plot(totalDPPK,totalDPPK*10,'k:','HandleVisibility','off')
plot(totalDPPK,totalDPPK*0.1,'k:','HandleVisibility','off')
plot([5*10^-5, 5*10^-5],[5*10^-8,1*10^0],'LineWidth',2)

legend('SDR m = 0','SDR m = 1','SDR m = 4','5*10^{-5} m/s line','Location','northwest')
ylim([5*10^-8,1*10^0])
xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')

%% Look at error factor

medianEF_m0 = median(NMR_K_m0_errorFactor)
medianEF_m1 = median(NMR_K_m1_errorFactor)
medianEF_m4 = median(NMR_K_m4_errorFactor)

            