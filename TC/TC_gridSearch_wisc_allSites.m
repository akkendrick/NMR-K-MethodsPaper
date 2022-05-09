% Run Timur-Coates Model Estimates
% Range over pairs of m and n values
 close all
 clear

load enso 

siteList = [{"Site1-WellG5"},{"Site1-WellG6"},{"Site2-WellPN1"},{"Site2-WellPN2"}];

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



for j = 1:length(siteList)
    currentFitMatrix = [];
    currentErrorMatrix = [];
    baseName = siteList{j};

    mTemp = zeros(1,length(cutoff));
    nTemp = zeros(1,length(cutoff));
    cTemp = zeros(1,length(cutoff));
    errorTemp = zeros(1,length(cutoff));
    meanFactorTemp = zeros(1,length(cutoff));
    indexQuotientTemp = zeros(1,length(cutoff));
    
    %Compute statistics for each cutoff, where the cutoffs are derived from
    %each T2linbins of the data
    
    %parfor i = 1:length(cutoff)
    for i = 1


        [K,z,T2dist,T2logbins,phi,kTC_best,bestFitMatrix,totalError,meanErrorFactor,...
            medianErrorFactor,indexQuotientLog] = computeTCperm(baseName,n,m,cutoff(i),figureson);

        mTemp(i) = bestFitMatrix(3);
        nTemp(i) = bestFitMatrix(2);
        cTemp(i) = bestFitMatrix(1);
        kTC{i} = kTC_best;
        
        %errorTemp(i) = log10(MAEError); %Why take log10 here?
        errorTemp(i) = totalError;
        meanFactorTemp(i) = meanErrorFactor;
        medianFactorTemp(i) = medianErrorFactor;

    end
    
    indexQuotient{j} = 10.^indexQuotientLog;
    
    DPP_K{j} = K;
    
    totalkTC{j} = kTC;
    totalmMatrix(j,:) = mTemp;
    totalnMatrix(j,:) = nTemp;
    totalcMatrix(j,:) = cTemp;
    totalErrorMatrix(j,:) = errorTemp;
    totalmeanErrorFactorMatrix(j,:) = meanFactorTemp;
    totalmedianErrorFactorMatrix(j,:) = medianFactorTemp;

    TC_bs{j} = TC_b(DPP_K{j},phi,indexQuotient{j});
end
toc



save('TC_33_out.mat','totalmMatrix','totalnMatrix','totalcMatrix',...
'totalErrorMatrix','totalmeanErrorFactorMatrix','totalmedianErrorFactorMatrix','cutoff','siteList','n','m')

%%
% Just run this part to save time
% Only works if using multiple cutoffs, otherwise there is only one data point
%load('optimalCutoffTable_n2_m1_RMSE_2000.mat')
%cutoff = (20:2:800)*10^-3;
%totalmeanErrorFactorMatrix = totalErrorFactorMatrix;
load('TC_33_out.mat')

smoothWindow = 5;

for kk = 1:length(siteList)
    
    dataSmooth = smooth(totalmedianErrorFactorMatrix(kk,:),smoothWindow);
    [minVal, index] = min(totalmedianErrorFactorMatrix(kk,:));

    minVals(kk) = minVal;
    minIndicies(kk) = index;
    
    cutoff_ms = cutoff*10^3;
    
    figure(1)
    subplot(2,2,kk)
    hold on
    plot(cutoff_ms,totalmedianErrorFactorMatrix(kk,:),'.')
 %   plot(cutoff_ms, dataSmooth,'LineWidth',2)
    plot(cutoff_ms(index),dataSmooth(index),'k*','MarkerSize',10)
    legend(siteList{kk})
    grid on
    box on
    xlabel('Cutoff (ms)')
    %ylabel('Mean Average Error (m/day)')
    ylabel('Median K Error Factor')
    set(gca,'FontSize',14)
    
  
    bestCutoff(kk,1) = cutoff_ms(index);
    bestC(kk,1) = totalcMatrix(kk,index);
    
    bestErrorMatrix(kk,1) = totalErrorMatrix(kk,index);
    bestMeanErrorFactorMatrix(kk,1) = totalmeanErrorFactorMatrix(kk,index);
    bestMedianErrorFactorMatrix(kk,1) = totalmedianErrorFactorMatrix(kk,index);
    
    ylim([0,10])
    xlim([0,max(cutoff_ms)])
            
%     edges = [0 0.1 0.2 0.5 1 2 3 4 5 10];
% 
%     figure(2)
%     subplot(2,2,kk)
%     histogram(totalErrorMatrix(kk,:),edges)
end
% 
% bestMeanErrorFactorMatrix
% bestMedianErrorFactorMatrix
% bestCutoff
% bestC
% 

