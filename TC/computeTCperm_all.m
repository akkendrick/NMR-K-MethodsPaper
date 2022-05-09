function [c_boot,kTC_best,medianErrorFactor] = computeTCperm_all(Kall,Zall,T2dist,T2logbins,logPhi,n,m,cutoff)

Dk = Kall;
z_dk = Zall;

origT2dist = T2dist;

zT2dist = T2dist(:,1);
T2dist = T2dist(:,2:end);

T2linbins = 10.^T2logbins;

logK = log(Kall);

% Match depths of permeability measurements from DPP data
% with the depths in the new processed data -- need to find nearest
% neighbors, exlude NMR points that do not have corresponding k
% measurement. 

errorz = zeros(1,length(Dk));
ind = zeros(1,length(Dk));

for i = 1:length(Dk)
   [errorz(i), ind(i)] = min(abs(z_dk(i) - zT2dist)); 
end

%% Run bootstrap
z = zT2dist(ind); 

% Estimate T-C parameters at DPP K intervals
% Specify cutoff between FFI and BVI in ms
[~, cutoffBin] = min(abs(T2linbins - cutoff));

filtT2dist = T2dist(ind,:);

% Signal amplitude is calibrated to sum to porosity
% sum(boundT2dist) + sum(freeT2dist) = phi

BVI = zeros(1,length(ind));
FFI = zeros(1,length(ind));
nBins = length(T2logbins);
for k = 1:length(ind)
    
    boundT2dist(1:cutoffBin) = filtT2dist(k,1:cutoffBin);
    boundT2dist(cutoffBin+1:nBins) = 0;

    freeT2dist(1:cutoffBin) = 0;
    freeT2dist(cutoffBin+1:nBins) = filtT2dist(k,cutoffBin+1:nBins);
    
    BVI(k) = sum(boundT2dist);
    
    FFI(k) = sum(freeT2dist);     
    
    if BVI(k) == 0
        BVI(k) = (randi(100)/1000)^randi(100);
    end
    
    if FFI(k) == 0
        FFI(k) = randi(100)*randi(100);
    end
    
end

lkTC = @(c,m,n,lPhi,logFrac) log10(c) + m*lPhi + n*(logFrac);

% Now that we have estimated BVI and FFI via cutoff, use bootstrap to
% estimate empirical parameters
indexQuotientLog = log10(FFI./BVI)';

Nboot = 2000;
%

[c_boot, n_boot, m_boot] = bootstrap_fun([indexQuotientLog, logPhi, logK], Nboot, n, m);

disp('Median prefactor for TC is:')
median_c = median(c_boot)

median_n = median(n_boot);
median_m = median(m_boot);

lkTC_best = lkTC(median_c,median_m,median_n,logPhi,indexQuotientLog);
kTC_best = 10.^lkTC_best;

disp('Median K diff factor for TC is:')
medianErrorFactor = median(estimateKdiffFactor(Kall,kTC_best,1))





end

