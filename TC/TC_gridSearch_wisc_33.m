% Run Timur-Coates Model Estimates
% Range over pairs of m and n values
close all
clear

load enso 

% siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'} {'Site1-WellG5'} {'Site1-WellG5above'}...
%     {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

siteList = [{"Site1-WellG5"},{"Site1-WellG6"},{"Site2-WellPN1"},{"Site2-WellPN2"}];
%siteList = [{'Site2-WellPN2'}];

offsets = [0.75,0.95,0.75,0.75];

% Matrix Structure... Repeat for each site
%       Bootstrap           
% c
% m
% n


cutoff = 33*10^-3;

m = 1;
n = 2;

figureson = 0;
wDirect = 1;
currentRow = -3;
%%
tic
totalErrorMatrix = zeros(length(siteList),length(cutoff));
totalmeanErrorFactorMatrix = zeros(length(siteList),length(cutoff));
totalmedianErrorFactorMatrix = zeros(length(siteList),length(cutoff));
totalmMatrix = zeros(length(siteList),length(cutoff));
totalnMatrix = zeros(length(siteList),length(cutoff));
totalcMatrix = zeros(length(siteList),length(cutoff));

for j = 1:length(siteList)
    currentFitMatrix = [];
    currentErrorMatrix = [];
    baseName = siteList{j};

    mTemp = zeros(1,length(cutoff));
    nTemp = zeros(1,length(cutoff));
    cTemp = zeros(1,length(cutoff));
    errorTemp = zeros(1,length(cutoff));
    meanFactorTemp = zeros(1,length(cutoff));
    
    %Compute statistics for each cutoff, where the cutoffs are derived from
    %each T2linbins of the data
    
    for i = 1:length(cutoff)


        [K,z,T2dist,T2logbins,kTC_best,bestFitMatrix,totalError,meanErrorFactor,medianErrorFactor,indexQuotient{j}] = computeTCperm(baseName,n,m,cutoff,figureson);

        mTemp(i) = bestFitMatrix(3);
        nTemp(i) = bestFitMatrix(2);
        cTemp(i) = bestFitMatrix(1);
        kTC{i} = kTC_best;
        
        %errorTemp(i) = log10(MAEError); %Why take log10 here?
        errorTemp(i) = totalError;
        meanFactorTemp(i) = meanErrorFactor;
        medianFactorTemp(i) = medianErrorFactor;

    end
    
    DPP_K{j} = K;
    
    totalkTC{j} = kTC;
    totalmMatrix(j,:) = mTemp;
    totalnMatrix(j,:) = nTemp;
    totalcMatrix(j,:) = cTemp;
    totalErrorMatrix(j,:) = errorTemp;
    totalmeanErrorFactorMatrix(j,:) = meanFactorTemp;
    totalmedianErrorFactorMatrix(j,:) = medianFactorTemp;

end
toc

save('TC_33_out.mat','totalmMatrix','totalnMatrix','totalcMatrix',...
'totalErrorMatrix','totalmeanErrorFactorMatrix','totalmedianErrorFactorMatrix','cutoff','siteList','n','m')


