function [Kfilt,z,b_boot,k_boot,totalErrorEstimate] = computeProfile_Kfilt(site,n,m,Kcutoff,figureson,wDirect,saveData)

[T2dist, T2logbins, nmrName] = loadRawNMRdata(site);

T2linbins = 10.^T2logbins;

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(nmrName); 


T2depths = T2dist(:,1);
T2data = T2dist(:,2:end);

for jj = 1:length(z)
    [minVal, index] = min(abs(T2depths - z(jj)));
    filtIndices(jj) = index;
    filtT2dist(jj,:) = T2data(index,:);
    filtT2depths(jj) = T2depths(index);
end

logSumEch = log10(SumEch);

bestFitMatrix = zeros(3,3);
totalErrorEstimate = zeros(1,3);

k_bootstrap = [];
k_mcmc = [];
k_direct = [];

[Ksort, ind] = sort(K);
T2ML_sort = T2ML(ind);
logT2ML_sort = logT2ML(ind);
logPhi_sort = logPhi(ind);
logK_sort = logK(ind);
phi_sort = phi(ind);

Kind = Ksort > Kcutoff;
Kfilt = Ksort(Ksort > Kcutoff);
logT2ML_filt = logT2ML_sort(Kind);
logPhi_filt = logPhi_sort(Kind);
logK_filt = logK_sort(Kind);
T2ML_filt = T2ML_sort(Kind);
phi_filt = phi_sort(Kind);


% Make sure we have valid data
if isempty(Kfilt)
    Kfilt = NaN;
    k_boot = NaN;
    z = NaN;
    b_boot = NaN;
    totalErrorEstimate = NaN; 
else
    %%
    % Use discretize to group data into bins
    % Assign each bin a weight based on the relative freq compared to the total
    % Use datasample with weights to try and obtain an even distribution, plot
    % hist for mult realizations, check results


    Nboot =  2000; % number of bootstrap samples

    % Takes [log10(T2ML), log10(K)] or [log10(T2ML), log10(phi), log10(K)] as a
    % single matrix


    [b_boot, n_boot, m_boot] = bootstrap_fixed_m_n([logT2ML_filt, logPhi_filt, logK_filt], Nboot, n, m);   % m, n fixed

    if figureson ==1 
        bs = log10(b_boot); 
        graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
    else
        bs = log10(b_boot); 
    end

    median_b = median(b_boot)
    std_b = std(b_boot);

    % Now compute k values as a function of depth
    k_boot = median_b*(phi_filt.^m).*(T2ML_filt).^n;

    bestFitMatrix(1,1) = median_b;
    bestFitMatrix(2,1) = m;
    bestFitMatrix(3,1) = n;

    totalErrorEstimate(1) = computeError(Kfilt, k_boot);
    totalErrorEstimate(2) = median(estimateKdiffFactor(Kfilt,k_boot,1))
    totalErrorEstimate(3) = std_b


    %graph_correlations([b_boot, m_boot], 2, {'log_{10}(b)', 'm'}, 1, 0)

    if saveData ==1 
        save(strcat(saveName,'_bootstrap_n_m_var.mat'),'bs','n_boot','m_boot')
    end
end

end

