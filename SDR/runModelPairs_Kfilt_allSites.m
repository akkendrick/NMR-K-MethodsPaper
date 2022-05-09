% Range over pairs of m and n values
close all
clear

%Wisconsin Data
siteList = [{'Site1-WellG5'}...
   {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
    
% Maurer and Knight
%  siteList = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%    'dpnmr_leque_east','dpnmr_leque_west'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fitting variables and K range
m = 1;
n = 2;
Kcutoff = 5*10^-5;
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
Kind = allK < Kcutoff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allK_small = allK(Kind);
allK_filt = allK(Kind);
allT2ML_filt = allT2ML(Kind);
allPhi_filt = allPhi(Kind);

alllogK_filt = alllogK(Kind);
alllogT2ML_filt = alllogT2ML(Kind);
alllogPhi_filt = alllogPhi(Kind);

    
[b_boot, n_boot, m_boot] = bootstrap_fixed_m_n([alllogT2ML_filt, alllogPhi_filt, alllogK_filt], Nboot, n, m);   % m, n fixed

median_b_small = median(b_boot)
NMR_K_small = median_b_small*(allPhi_filt.^m).*(allT2ML_filt).^n;

%% Now do for opposite direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fitting variables and K range
m = 1;
n = 2;
Kcutoff = 5*10^-5;
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

allK_large = allK(Kind);
allK_filt = allK(Kind);
allT2ML_filt = allT2ML(Kind);
allPhi_filt = allPhi(Kind);

alllogK_filt = alllogK(Kind);
alllogT2ML_filt = alllogT2ML(Kind);
alllogPhi_filt = alllogPhi(Kind);

    
[b_boot, n_boot, m_boot] = bootstrap_fixed_m_n([alllogT2ML_filt, alllogPhi_filt, alllogK_filt], Nboot, n, m);   % m, n fixed

median_b_large = median(b_boot)
NMR_K_large = median_b_large*(allPhi_filt.^m).*(allT2ML_filt).^n;

%% Now do for all the points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

median_b_all = median(b_boot)
NMR_K_all = median_b_all*(allPhi_filt.^m).*(allT2ML_filt).^n;
%% Now plot results
totalDPPK = allK;

 figure()

hold on
    
grid on 
box on


scatter(allK_large, NMR_K_large,60,'Filled')
scatter(allK_small, NMR_K_small,60,'Filled')
scatter(totalDPPK, NMR_K_all,60,'Filled')

plot(totalDPPK,totalDPPK,'k','LineWidth',2,'HandleVisibility','off')
plot(totalDPPK,totalDPPK*10,'k:','HandleVisibility','off')
plot(totalDPPK,totalDPPK*0.1,'k:','HandleVisibility','off')
plot([5*10^-5, 5*10^-5],[5*10^-8,1*10^0],'LineWidth',2)

legend('SDR Large K','SDR Small K','SDR All K','5*10^{-5} m/s line','Location','northwest')
ylim([5*10^-8,1*10^0])
xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')


            