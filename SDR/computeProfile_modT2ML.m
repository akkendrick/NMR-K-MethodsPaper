function [K,z,k_bootstrap,k_mcmc,k_direct,bestFitMatrix,b_boot,T2ML,phi,totalErrorEstimate,finalDepths] = computeProfile_modT2ML(site,n,m,figureson,wDirect,saveData)

[T2dist, T2logbins, nmrName] = loadRawNMRdata(site);

T2linbins = 10.^T2logbins;

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(nmrName); 

depthsAll = z;
finalDepths = depthsAll;

% Filter data by depth
if (site == "Site1-WellG6")


    depthCutoff = 5.8;

    K = K(depthsAll>depthCutoff);
    phi = phi(depthsAll>depthCutoff);
    T2ML = T2ML(depthsAll>depthCutoff);
    SumEch = SumEch(depthsAll>depthCutoff);
    logK = logK(depthsAll>depthCutoff);
    logT2ML = logT2ML(depthsAll>depthCutoff);
    logPhi = logPhi(depthsAll>depthCutoff);
    T2dist = T2dist(T2dist(:,1)>depthCutoff,:);
    finalDepths = depthsAll(depthsAll>depthCutoff);
    
elseif (site == "Site1-WellG5")


    depthCutoff = 4;
    
    K = K(depthsAll>depthCutoff);
    phi = phi(depthsAll>depthCutoff);
    T2ML = T2ML(depthsAll>depthCutoff);
    SumEch = SumEch(depthsAll>depthCutoff);
    logK = logK(depthsAll>depthCutoff);
    logT2ML = logT2ML(depthsAll>depthCutoff);
    logPhi = logPhi(depthsAll>depthCutoff);
    T2dist = T2dist(T2dist(:,1)>depthCutoff,:);  
    finalDepths = depthsAll(depthsAll>depthCutoff);


end


T2depths = T2dist(:,1);
T2data = T2dist(:,2:end);
logSumEch = log10(SumEch);

for jj = 1:length(z)
    [minVal, index] = min(abs(T2depths - z(jj)));
    filtIndices(jj) = index;
    filtT2dist(jj,:) = T2data(index,:);
    filtT2depths(jj) = T2depths(index);
end



bestFitMatrix = zeros(3,3);
totalErrorEstimate = zeros(1,3);

k_bootstrap = [];
k_mcmc = [];
k_direct = [];
%% Bootstrap

% Provide bootstrap an evenly weighted distribution points
% Idea:

%edges = [10^-6,5*10^-6,10^-5,5*10^-5,10^-4,5*10^-4,10^-3];
%edges = [10^-6,10^-4,10^-3];


% edges = [10^-6,5*10^-5,10^-3];
% discretizedK = discretize(K,edges);
% 
% nPts = length(K);
% %weights(discretizedK == 1) =  1-sum(discretizedK == 1)/nPts;
% % weights(discretizedK == 1) =  1;
% % 
% % weights(discretizedK == 2) =  1-sum(discretizedK == 2)/[[nPts;
% 
% [targetPoints, index] = min([sum(discretizedK == 1), sum(discretizedK == 2)]);
% 
% sampledK_2 = datasample(K(discretizedK == 2),targetPoints,'Replace',false); 
% sampledK_1 = datasample(K(discretizedK == 1),targetPoints,'Replace',false);
% 
% sampledK = [sampledK_2; sampledK_1];
% 
% discretizedSample = discretize(sampledK,edges);
% hist(discretizedSample)
% 
% % Now sampledK has the small K and a random selection of high K samples
% for kk = 1:length(sampledK)
%     % Take the min just in case we have multiple K values that are the same
%     filtIndex = min(find(K == sampledK(kk)));
%     logT2ML_filt(kk) = logT2ML(filtIndex);
%     logPhi_filt(kk) = logPhi(filtIndex);
%     logK_filt(kk) = logK(filtIndex);
% end
% 
% logT2ML_filt = logT2ML_filt';
% logPhi_filt = logPhi_filt';
% logK_filt = logK_filt';

%% %Compute different T2ML metric
for jj = 1:length(T2ML)
    sliceT2 = filtT2dist(jj,:);
    [pks{jj},locs{jj},widths{jj},proms{jj}] = findpeaks(sliceT2);

    if length(pks{jj}) > 1
        T2peakAvg(jj) =  mean(T2linbins(locs{jj}));
        currentProm = proms{jj};
        currentWidth = widths{jj};
        currentLocs = locs{jj};

        promRatio(jj,1) = currentProm(1)/currentProm(2);
        widthRatio(jj,1) = currentWidth(1)/currentWidth(2);

        T2peakMax(jj) = T2linbins(currentLocs(2));
        T2peakMin(jj) = T2linbins(currentLocs(1));
    else
        
        T2peakAvg(jj) = T2linbins(locs{jj});
        promRatio(jj,1) = NaN;
        widthRatio(jj,1) = NaN;
        T2peakMax(jj) = T2linbins(locs{jj});
        T2peakMin(jj) = NaN;
    end

    % Compute summed portions of sliceT2 based on transition between large
    % and small pores
    minLoc = islocalmin(sliceT2);
    [maxVal, maxInd] = max(sliceT2);
    
    stepLocalMax = islocalmax(sliceT2);
    T2localMax = T2linbins(stepLocalMax);
    
    diffPeaks = T2peakMax - T2peakMin;
    
    minInd = find(minLoc == 1);

    % Take minVals less than max
    goodMinLoc = minInd(minInd < maxInd);
    
    % Take the local min with the largest relax time
    goodMinLoc = max(goodMinLoc);
    allGoodMinLoc{jj,1} = goodMinLoc;

    T2hm_quant = 0;
    sliceT2_filt = ones(1,length(sliceT2));
    for kk = 1:1:length(sliceT2)
        %Assuming T2bins in linear space
         T2hm_quant = T2hm_quant + sliceT2(kk)/T2linbins(kk);
         %Is sliceT2 really the quantitfy we want here?
         if sliceT2(1,kk) ~= 0
             sliceT2_filt(1,kk) = sliceT2(1,kk);
         end
    end
    
    T2gm_quant = cumprod(sliceT2_filt,2);
    
    T2am(jj,1) = sum(sliceT2)/length(sliceT2);
    T2gm(jj,1) = T2gm_quant(end)^(1/length(sliceT2));
    T2hm(jj,1) = sum(sliceT2)/T2hm_quant;
    T2lm(jj,1) = 10.^(sum(T2logbins.*sliceT2)./phi(jj));
    
    if ~isempty(goodMinLoc)
        shortSignal(1:goodMinLoc) = sliceT2(1:goodMinLoc); 
        shortSignal(goodMinLoc+1:100) = 0;
        
        [maxVal, indShortMax(jj,1)] = max(shortSignal);
        shortT2Max(jj,1) = T2linbins(indShortMax(jj,1));
        
        longSignal(1:goodMinLoc) = 0; 
        longSignal(goodMinLoc+1:100) = sliceT2(goodMinLoc+1:end);
        
        [maxVal, indLongMax(jj,1)] = max(longSignal);
        longT2Max(jj,1) = T2linbins(indLongMax(jj,1));
        
        longSignals(1:100,jj) = longSignal;
        shortSignals(1:100,jj) = shortSignal;
%         for j = 2:length(T2linbins)
%         
%             avgShort = (shortSignal(j)+shortSignal(j-1))/2; 
%             avgLong = (longSignal(j)+longSignal(j-1))/2; 
% 
%             shortTemp(j) = avgShort*T2spacing(j-1);
%             longTemp(j) = avgLong*T2spacing(j-1);  
%         end
    
        shortSummedSignal(jj) = sum(shortSignal);
        longSummedSignal(jj) = sum(longSignal);
        
        allMinLoc(jj,1) = goodMinLoc;
        allT2min(jj,1) = T2linbins(goodMinLoc);
                
        T2lm_short(jj,1) = 10.^(sum(T2logbins(1:goodMinLoc).*sliceT2(1:goodMinLoc))./shortSummedSignal(jj));
        T2lm_long(jj,1) = 10.^(sum(T2logbins(goodMinLoc+1:end).*sliceT2(goodMinLoc+1:end))./longSummedSignal(jj));
            
%         modT2ML(jj,1) = (longSummedSignal(jj)/phi(jj,1))*T2lm_long(jj,1) + (shortSummedSignal(jj)/phi(jj,1))*T2lm_short(jj,1);
%         modT2ML(jj,1) = (longSummedSignal(jj) + shortSummedSignal(jj))/(longSummedSignal(jj)/T2lm_long(jj,1)+...
%             shortSummedSignal(jj)/T2lm_short(jj,1));
%           modT2ML(jj,1) = ((T2lm_long(jj,1)^longSummedSignal(jj))*...
%             T2lm_short(jj,1)^shortSummedSignal(jj))^(1/phi(jj));
    else
      
        longSignal(1:100) = sliceT2(1:end);
        longSignals(1:100,jj) = longSignal;
        shortSignals(1:100,jj) = zeros(1,100);
        
%         for j = 2:length(T2linbins)
%             avgLong = (longSignal(j)+longSignal(j-1))/2; 
%             longTemp(j) = avgLong*T2spacing(j-1);  
%         end
        
        shortSummedSignal(jj) = 0;
        longSummedSignal(jj) = sum(longSignal);
        allMinLoc(jj,1) = NaN;
        T2lm_short(jj,1) = NaN;
        allT2min(jj,1) = NaN;
        T2lm_long(jj,1) = T2ML(jj);
        
        T2_bimodal(jj) = 0;
        
        [maxVal, indLongMax(jj)] = max(longSignal);
        longT2Max(jj,1) = T2linbins(indLongMax(jj));
        shortT2Max(jj,1) = 0;
        indShortMax(jj,1) = 1;
    
%         modT2ML(jj,1) = ((T2lm_long(jj,1)^longSummedSignal(jj)))^(1/phi(jj));

   end

   
    
end

%modT2ML = T2am;
%modT2ML = T2hm;
%modT2ML = T2gm;
%modT2ML = T2lm;
modT2ML = T2ML;

%Fix short/long summed signal! Bin spacing is nonlinear!
T2spacing = diff(T2linbins);
largeModeSum(jj,1) = longSummedSignal(jj);
smallModeSum(jj,1) = shortSummedSignal(jj);
allT2cumsum(jj,:) = cumsum(sliceT2,'reverse');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logT2ML_filt = log10(modT2ML);
%logT2ML_filt = logT2ML;
logPhi_filt = logPhi;
logK_filt = logK;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edges for Wisconsin
% Only look at "small values of K"
%lowestKEdge = 6*10^-5;
%midKEdge = 10^-4; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edges for Maurer and Knight
% Only look at "small values of K"
% lowestKEdge = 10^-5;
% midKEdge = 5*10^-4; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %filtInd = find(K < lowestKEdge);
% %filtInd = find(K < midKEdge & K > lowestKEdge);
% filtInd = find(K > midKEdge);
% 
% logT2ML_filt = logT2ML(filtInd);
% logPhi_filt = logPhi(filtInd);
% logK_filt = logK(filtInd);

%%
% Use discretize to group data into bins
% Assign each bin a weight based on the relative freq compared to the total
% Use datasample with weights to try and obtain an even distribution, plot
% hist for mult realizations, check results


Nboot =  2000; % number of bootstrap samples

% Takes [log10(T2ML), log10(K)] or [log10(T2ML), log10(phi), log10(K)] as a
% single matrix
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt,lp, kk], Nboot);         % m, n can vary
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt, lp, kk], Nboot, n);        % m can vary
%  [b_boot, n_boot, m_boot] = bootstrap_fun_mb([logT2ML, logK], Nboot);    % n can vary
  
%  [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot);    % n and m can vary

% if isempty(n) && isempty(m)
%     [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot);
% elseif ~isempty(n) && isempty(m)
%     [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot, n);
% else
%     [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot, n, m);   % m, n fixed
% end

if isempty(n) && isempty(m)
    [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML_filt, logPhi_filt, logK_filt], Nboot);
elseif ~isempty(n) && isempty(m)
    [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML_filt, logPhi_filt, logK_filt], Nboot, n);
else
    [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML_filt, logPhi_filt, logK_filt], Nboot, n, m);   % m, n fixed
end

if figureson ==1 
    bs = log10(b_boot); 
    graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
else
    bs = log10(b_boot); 
end


meanb = mean(b_boot);
sortb = sort(b_boot); 
% blo = sortb(50);
% bhi = sortb(1950);

median_b = median(b_boot)
std_b = std(b_boot);

% BIG Q: Use mean or median? Maybe should be using median here... 

% Now compute k values as a function of depth
if isempty(m_boot) && isempty(n_boot)
    k_bootstrap = median_b*(phi.^m).*(T2ML).^n;
    
    bestFitMatrix(1,1) = median_b;
    bestFitMatrix(2,1) = m;
    bestFitMatrix(3,1) = n;
    
    totalErrorEstimate(1) = computeError(K, k_bootstrap);
else
    if isempty(n_boot) && ~isempty(m_boot)
        median_n = median(n_boot);
        k_bootstrap = median_b*(phi.^m).*(T2ML).^median_n;
       
        bestFitMatrix(1,1) = median_b;
        bestFitMatrix(2,1) = m;
        bestFitMatrix(3,1) = median_n;
        
        totalErrorEstimate(1) = computeError(K, k_bootstrap);

    else
        median_m = median(m_boot);
        median_n = median(n_boot);
        k_bootstrap = median_b*(phi.^median_m).*(T2ML).^median_n;
        
        bestFitMatrix(1,1) = median_b;
        bestFitMatrix(2,1) = median_m;
        bestFitMatrix(3,1) = median_n;
        
        totalErrorEstimate(1) = computeError(K, k_bootstrap);
        totalErrorEstimate(2) = median(estimateKdiffFactor(K,k_bootstrap,1))
        totalErrorEstimate(3) = std_b
    end
end

%graph_correlations([b_boot, m_boot], 2, {'log_{10}(b)', 'm'}, 1, 0)

if saveData ==1 
    save(strcat(saveName,'_bootstrap_n_m_var.mat'),'bs','n_boot','m_boot')
end
%% Basic solving for b for fixed n, m
    % Given m and n, we can solve directly for b -> b = log(k) - m*log(phi) -
    % n*log(T2ML). This is the 'direct' method.

if wDirect == 1
    mDirect = bestFitMatrix(2,1);
    nDirect = bestFitMatrix(3,1);
    
    C = @(m, n, lt, lp) m*lp + n*lt; 

    bdir_n2 = logK - C(mDirect, nDirect, logT2ML, logPhi); 

    logb_mean = mean(bdir_n2);
    b_mean = 10.^logb_mean;
    bs_basic = 10.^(bdir_n2);


    logb_median = median(bdir_n2);
    median_b = 10.^logb_median;

    k_direct = b_mean*(phi.^mDirect).*(T2ML).^nDirect;

    bestFitMatrix(1,2) = median_b;
    bestFitMatrix(2,2) = mDirect;
    bestFitMatrix(3,2) = nDirect;

    totalErrorEstimate(2) = computeError(K, k_direct);
end

if saveData == 1
    save(strcat(saveName,'_basic_solving_n_m_var.mat'),'bs_basic','nDirect','mDirect')
end

% %% MCMC for solution to various parameters
% % Markov Chain Monte Carlo using Metropolis-Hastings algorithm. Assumes
% % Bayes' theorem. Parameters Niter (number of iterations) and stepsize may
% % need to be adjusted for good covergence, although with this data set the
% % values shown perform reasonably well.
% Niter= 1e6; 
% stepsize = 0.8; 
% 
% %%%%%%%%%%%%%%%%%%%% MCMC hps gives summary parameters (average, CI, etc.) Pick
% % either this one or simpler MCMC algorithm to compare with other results.
% % NOTE: if the first MCMC algorthm is used, remember when comparing results
% % between different methods that the others use fixed values for m and n,
% % while this allows both to vary. This will inevitably affect the range in
% % b obtained. 
% 
% if isempty(m) || isempty(n)
% 
%     %%%%%%%%%%%%%%%%%%% MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
%     [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_full(K, T2ML, phi, ...
%         z, Niter, stepsize, figureson);
%     T2B_mcmc  = paramhats(1,:);
%     blog_mcmc = paramhats(2,:); 
%     n_mcmc = paramhats(3,:);
%     m_mcmc = paramhats(4,:);
%     sig_mcmc = paramhats(5,:);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Compute b statistics from variable n and m data
%     b_mcmc = 10.^blog_mcmc;
%     b_mean = mean(b_mcmc);
%     b_median = median(b_mcmc);
%     
%     n_mean = mean(n_mcmc);
%     n_median = median(n_mcmc);
%     
%     m_mean = mean(m_mcmc);
%     m_median = median(m_mcmc);
%     
%     k_mcmc = b_median*(phi.^m_median).*(T2ML).^n_median;
%      
%     bestFitMatrix(1,3) = b_median;
%     bestFitMatrix(2,3) = m_median;
%     bestFitMatrix(3,3) = n_median;
% 
%     totalErrorEstimate(3) = computeError(K, k_mcmc);
% 
% else 
% 
%     %%%%%%%%%%%%%% MCMC for b and data error only. n and m is fixed 
% %     [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_full(K, T2ML, phi, ...
% %         z, Niter, stepsize, figureson);
% %     T2B_mcmc  = paramhats(1,:);
% %     blog_mcmc = paramhats(2,:); 
% %     n_mcmc = paramhats(3,:);
% %     m_mcmc = paramhats(4,:);
% %     sig_mcmc = paramhats(5,:);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     
%     
%     [blog_mcmc, sig_mcmc, likes, accept_rat] = mcmc_nmr_bsig (K, T2ML, phi, z, m, n, ...
%        Niter, stepsize, figureson); 
% 
%     
%  %   Compute b statistics from fixed n and m data
%     b_mcmc = 10.^blog_mcmc;
%     b_mean = mean(b_mcmc);
%     b_median = median(b_mcmc);
%     
%     k_mcmc = b_median*(phi.^m).*(T2ML).^n;
%     
%     bestFitMatrix(1,3) = b_median;
%     bestFitMatrix(2,3) = m;
%     bestFitMatrix(3,3) = n;
%     
%     totalErrorEstimate(3) = computeError(K, k_mcmc);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     params = [blog_mcmc; sig_mcmc];
%    %  graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
%     graph_correlations(params,2,{'MCMC log_{10}(b)','MCMC \sigma'},1,0)
%    
%     
%     
% end
% 
% if saveData == 1
%     save(strcat(saveName,'_MCMC_n_m_var.mat'),'blog_mcmc','n_mcmc','m_mcmc')
% end
% 
% % 
% % if figureson == 1
% %     graph_correlations(paramhats, 2, {'T_B', 'log_{10}(b)', 'n', 'm', '\sigma'}, 0, 0)
% %     all_lKpreds = zeros(length(T2ML), length(blog_mcmc));   
% % %     for k = 1:length(blog_mcmc)
% % %         bbm = [blog_mcmc(k), sig_mcmc(k)];
% % %         [~,all_lKpreds(:,k)] = NMRfun2(bbm, K, phi, T2ML, m, n);
% % %     end
% % %     figure;
% % %     hold on
% % %     for k = 1:size(all_lKpreds,1)
% % %         dpk = all_lKpreds(k,:);
% % %         [h,x] = hist(dpk, 40);
% % %         t2m = repmat(log10(T2ML(k)), 1, length(x));
% % %         scatter(t2m, x, [], h);
% % %     end
% % %     scatter(log10(T2ML), logK, 'ok', 'filled')
% % end
% 
% %m_mcmc = 0;


end

