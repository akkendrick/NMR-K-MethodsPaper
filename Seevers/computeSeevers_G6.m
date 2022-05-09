function [b_boot,k_boot,totalErrorEstimate] = computeSeevers_G6(Kall,logK,...
    logT2ML,logPhi,n,m,depthsAll,depthCutoff)

% Compute Seevers Models for G6

% Make T2B vector
T2B = ones(length(depthsAll),1);
T2B(depthsAll > depthCutoff) = 3.048;
T2B(depthsAll < depthCutoff) = 2.433;

% interpT2B = 2.0680; % Fix T2B for other sites where we don't have data 
% % this is T2B from KGM model equation for T2B at 7 deg C estimated from
% % map of groundwater temperature
% % (https://pubs.usgs.gov/wsp/0520f/report.pdf) 
% % "Temperature of water available for industrial use in the United States: 
% % Chapter F in Contributions to the hydrology of the United States,
% % 1923-1924" by W.D. Collins
% % Equation for T2B from Dlugosch et al. 2013

phi = 10.^logPhi;
T2ML = 10.^logT2ML;
    
seeversT2 = (T2ML.^(-1) - T2B.^(-1)).^(-1);
logSeeversT2 = log10(seeversT2);

Nboot = 2000;


[b_boot, n_boot, m_boot] = bootstrap_fun([logSeeversT2, logPhi, logK], Nboot, n, m);   % m, n fixed

    
meanb = mean(b_boot);
disp('Median b is:')
median_b = median(b_boot)

% Now compute k values as a function of depth
median_m = median(m_boot);
median_n = median(n_boot);

k_boot = median_b*(phi.^m).*((seeversT2.^(-1)-T2B.^(-1)).^(-1)).^n;


disp('Median K diff factor is:')
totalErrorEstimate = median(estimateKdiffFactor(Kall, k_boot, 1))

end