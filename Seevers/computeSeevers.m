function [DK,z,T2dist,T2logbins,k_bootstrap,k_mcmc,k_direct,bestFitMatrix,...
    b_boot,idealbs,totalErrorEstimate,T2ML,T2B,phi] = computeSeevers(site,n,m,figureson,wDirect)

% Compute Seevers Models
%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
% baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Kansas_Wash_Data/';
Seevers_b = @(m,n,phi,T2ML,T2B,K) K ./ ((phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n);

[T2dist, T2logbins,nmrName] = loadRawNMRdata(site);

[d, DK, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(nmrName); 

phi_copy = phi;
T2ML_copy = T2ML;
depthsAll = z;

KcutoffLower = 0;
KcutoffHigher = 5*10^(5);

logSumEch = log10(SumEch);

bestFitMatrix = zeros(3,3);
totalErrorEstimate = zeros(1,3);

bestFitMatrix(:,2) = NaN;
totalErrorEstimate(2) = NaN;

k_bootstrap = [];
k_mcmc = [];
k_direct = [];

if (site == "Site1-WellG6")


    % Using T2B from Keating and Knight 2007
    % for quartz sand and ferrihydrite-coated sand

    depthCutoff = 5.8;
    T2B = 2.0680;

    DK = DK(depthsAll>depthCutoff);
    phi = phi(depthsAll>depthCutoff);
    T2ML = T2ML(depthsAll>depthCutoff);
    SumEch = SumEch(depthsAll>depthCutoff);
    logK = logK(depthsAll>depthCutoff);
    logT2ML = logT2ML(depthsAll>depthCutoff);
    logPhi = logPhi(depthsAll>depthCutoff);
    
%     T2B = ones(length(depthsAll),1);
%     T2B(depthsAll > depthCutoff) = 3.048;
%     T2B(depthsAll <= depthCutoff) = 2.433;
elseif (site == "Site1-WellG5")

    % Using T2B from Keating and Knight 2007
    % for quartz sand and ferrihydrite-coated sand

    depthCutoff = 4;
    T2B = 2.0680;

    DK = DK(depthsAll>depthCutoff);
    phi = phi(depthsAll>depthCutoff);
    T2ML = T2ML(depthsAll>depthCutoff);
    SumEch = SumEch(depthsAll>depthCutoff);
    logK = logK(depthsAll>depthCutoff);
    logT2ML = logT2ML(depthsAll>depthCutoff);
    logPhi = logPhi(depthsAll>depthCutoff);
    
%     T2B = ones(length(depthsAll),1);
%     T2B(depthsAll > depthCutoff) = 3.048;
%     T2B(depthsAll <= depthCutoff) = 2.433;
else
    %T2B = 3.048;
        
    T2B = 2.0680; % Fix T2B for other sites where we don't have data
    % this is T2B from KGM model equation for T2B at 7 deg C estimated from
    % map of groundwater temperature
    % (https://pubs.usgs.gov/wsp/0520f/report.pdf) 
    % "Temperature of water available for industrial use in the United States: 
    % Chapter F in Contributions to the hydrology of the United States,
    % 1923-1924" by W.D. Collins
    % Equation for T2B from Dlugosch et al. 2013
end

%% Plot T2B data for reference
% figure
% hold on
% plot(interpT2B, nmrDepths, '*')
% plot(T2B_peak, T2B_depth, '+')
% set(gca,'YDir','reverse')
% legend('Interpolated','Measured')
% grid on
% box on
% xlabel('T2B (s)')
% ylabel('Depth (m)')
% set(gca,'FontSize',14)

%% Implement Bootstrap
%seeversT2 = (T2ML.^-1 - interpT2B.^-1).^(-1);

seeversT2 = (T2ML.^(-1) - T2B.^(-1)).^(-1);
logSeeversT2 = log10(seeversT2);

Nboot = 2000;

if isempty(n) && isempty(m)
    [b_boot, n_boot, m_boot] = bootstrap_fun([logSeeversT2, logPhi, logK], Nboot);
else   
    [b_boot, n_boot, m_boot] = bootstrap_fun([logSeeversT2, logPhi, logK], Nboot, n, m);   % m, n fixed
end

if figureson ==1 
    bs = log10(b_boot); 
%    graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
else
    bs = log10(b_boot); 
end

meanb = mean(b_boot);

median_b = median(b_boot);

% Now compute k values as a function of depth
if isempty(m_boot) && isempty(n_boot)
    k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-T2Bavg.^(-1)).^(-1)).^n;

    %k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-interpT2B.^(-1)).^(-1)).^n;
    
    bestFitMatrix(1,1) = median_b;
    bestFitMatrix(2,1) = m;
    bestFitMatrix(3,1) = n;
    
    totalErrorEstimate(1) = computeError(DK, k_bootstrap);
    totalErrorEstimate(2) = mean(estimateKdiffFactor(DK, k_bootstrap, 1));

else
    if isempty(m_boot)
        median_n = median(n_boot);
        k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-T2Bavg.^(-1)).^(-1)).^n;
        %k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-interpT2B.^(-1)).^(-1)).^n;

        bestFitMatrix(1,1) = median_b;
        bestFitMatrix(2,1) = m;
        bestFitMatrix(3,1) = median_n;
        
        totalErrorEstimate(1) = computeError(DK, k_bootstrap);
        totalErrorEstimate(2) = mean(estimateKdiffFactor(DK, k_bootstrap, 1));

    else
        median_m = median(m_boot);
        median_n = median(n_boot);
        %k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-T2Bavg.^(-1)).^(-1)).^n;
        k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-T2B.^(-1)).^(-1)).^n;

        
        bestFitMatrix(1,1) = median_b;
        bestFitMatrix(2,1) = median_m;
        bestFitMatrix(3,1) = median_n;
        
        totalErrorEstimate(1) = computeError(DK, k_bootstrap);
        totalErrorEstimate(2) = mean(estimateKdiffFactor(DK, k_bootstrap, 1));
        
    end
end

idealbs = Seevers_b(m,n,phi,T2ML,T2B,DK);
%%
%  %%%%%%%%%%%%%%%%%%% MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
% Niter= 1e6; 
% stepsize = 0.8; 
%  
% [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_T2B(K, T2ML, interpT2B, phi, ...
%     z, Niter, stepsize, figureson);
% blog_mcmc = paramhats(1,:); 
% n_mcmc = paramhats(2,:);
% m_mcmc = paramhats(3,:);
% sig_mcmc = paramhats(4,:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Compute b statistics from variable n and m data
% b_mcmc = 10.^blog_mcmc;
% b_mean = mean(b_mcmc);
% b_median = median(b_mcmc);
% 
% n_mean = mean(n_mcmc);
% n_median = median(n_mcmc);
% 
% m_mean = mean(m_mcmc);
% m_median = median(m_mcmc);
% 
% k_mcmc = b_median*(phi.^m_median).*(T2ML).^n_median;
% 
% bestFitMatrix(1,3) = b_median;
% bestFitMatrix(2,3) = m_median;
% bestFitMatrix(3,3) = n_median;
% 
% totalErrorEstimate(3) = computeError(K, k_mcmc);




end