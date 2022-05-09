function [K,b_boot,k_boot,totalErrorEstimate]  = computeProfile_all(K,logK,logT2ML,logPhi,n,m,saveData)
%COMPUTEPROFILE_ALL Compute bootstrap of K data in a supplied vector
    
Nboot =  2000; % number of bootstrap samples

[b_boot, n_boot, m_boot] = bootstrap_fixed_m_n([logT2ML, logPhi, logK], Nboot, n, m);   % m, n fixed


median_b = median(b_boot)
std_b = std(b_boot);

phi = 10.^logPhi;
T2ML = 10.^logT2ML;

k_boot = median_b*(phi.^m).*(T2ML).^n;
totalErrorEstimate = median(estimateKdiffFactor(K,k_boot,1))

end

