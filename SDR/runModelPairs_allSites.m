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

m = 1;
n = 2;


saveData = 0;


for i = 1:length(siteList)
    siteName = siteList{i};
    [T2dist{i}, T2logbins{i}, nmrName{i}] = loadRawNMRdata(siteName);
    [d{i}, K{i}, T2ML{i}, phi{i}, z{i}, SumEch{i}, logK{i}, logT2ML{i}, logPhi{i}, SumEch_3s{i}, SumEch_twm{i}, ...
    SumEch_twm_3s{i}] = loadnmrdata2(nmrName{i}); 
    
    depthsAll = z{i}
    % Filter data by depth
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
    end
end

Kall = vertcat(K{:});
logK_all = vertcat(logK{:});
logT2ML_all = vertcat(logT2ML{:});
logPhi_all = vertcat(logPhi{:});
T2ML_all = vertcat(T2ML{:});
phi_all = vertcat(phi{:});
[K,b_boot,k_boot,totalErrorEstimate] = computeProfile_all(Kall,logK_all,...
    logT2ML_all,logPhi_all,n,m,saveData);
%%
disp('Median prefactor for SDR is:')
medianb = median(b_boot)

SDR_K = medianb.*(phi_all.^m).*(T2ML_all.^n);

disp('Median K diff factor for SDR is:')
errorEstimate = median(estimateKdiffFactor(Kall, SDR_K,1))
% save('SDR_T2lm_bestFit_092821_m1_n2.mat','siteList','m','n','matrixKey','totalbMatrix','totalmMatrix','totalErrorMatrix','totalKMatrix')
            
            