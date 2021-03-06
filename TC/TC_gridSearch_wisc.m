% Run Timur-Coates Model Estimates
% Range over pairs of m and n values
%  close all
%  clear

load enso 

% siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'} {'Site1-WellG5'} {'Site1-WellG5above'}...
%     {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

siteList = [{"Site1-WellG5"},{"Site1-WellG6"},{"Site2-WellPN1"},{"Site2-WellPN2"}];
%siteList = [{'Site2-WellPN2'}];

TC_b = @(K,phi,frac) K ./ ((phi).*(frac).^2);
offsets = [0.75,0.95,0.75,0.75];

% Matrix Structure... Repeat for each site
%       Bootstrap           
% c
% m
% n

%cutoff = [20 33 40 50 60 80 100 120 140 160 200 240 280 320 380 440 500 560 620]*10^-3;
%  cutoff_1 = (10:2:800)*10^-3;
%  cutoff_2 = (800:20:1200)*10^-3;
%  cutoff = [cutoff_1 cutoff_2];

% Using T2linbins from Wisc data
%cutoff = [0.00112174926178712,0.00122548246058477,0.00133880833476839,0.00146261396216966,0.00159786845269657,0.00174563053420842,0.00190705683989084,0.00208341096200154,0.00227607334285676,0.00248655208048162,0.00271649473350762,0.00296770121772206,0.00324213789521941,0.00354195296643983,0.00386949328468776,0.00422732273001957,0.00461824227333737,0.00504531190480970,0.00551187458565686,0.00602158241575669,0.00657842522108640,0.00718676178477373,0.00785135396622390,0.00857740397539565,0.00937059509399720,0.0102371361623543,0.0111838101801791,0.0122180274016706,0.0133478833405578,0.0145822221391302,0.0159307057936191,0.0174038897995025,0.0190133057553889,0.0207715516423368,0.0226923904333668,0.0247908578255075,0.0270833799343078,0.0295879018720907,0.0323240282164155,0.0353131764682914,0.0385787447013611,0.0421462947143596,0.0460437521205034,0.0503016249400461,0.0549532424070700,0.0600350158598163,0.0655867237416149,0.0716518230326417,0.0782777893301519,0.0855164885284808,0.0934245827944563,0.102063974100286,0.111502288771921,0.121813406845797,0.133078040377556,0.145384365229593,0.158828711282831,0.173516316475480,0.189562150571158,0.207091815104548,0.226242526549083,0.247164190402560,0.270020574536119,0.294990591358937,0.322269697928723,0.352071426158467,0.384629054217416,0.420197432556189,0.459054978792650,0.501505857073845,0.547882358973335,0.598547504560823,0.653897884004373,0.714366761948377,0.780427468967305,0.852597106642439,0.931440595263711,1.01757509560643,1.11167484078113,1.21447641256319,1.32678450799154,1.44947823805903,1.58351800910544,1.72995304056379,1.88992957790645,2.06469986507934]

% Estimating ideal cutoff from data to separate T2 modes (this is in s)
%cutoff = [0.0102371361623543,0.0111838101801791,0.0122180274016706,0.0133478833405578,0.0145822221391302,0.0159307057936191,0.0174038897995025,0.0190133057553889,0.0207715516423368,0.0226923904333668,0.0247908578255075,0.0270833799343078,0.0295879018720907,0.0323240282164155,0.0353131764682914,0.0385787447013611,0.0421462947143596,0.0460437521205034,0.0503016249400461,0.0549532424070700,0.0600350158598163,0.0655867237416149,0.0716518230326417,0.0782777893301519,0.0855164885284808,0.0934245827944563,0.102063974100286];


%cutoff = 42.1*10^-3; %using median of the min T2 time for bimodal distributions
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

