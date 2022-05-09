

sites = ["Site1-WellG5","Site1-WellG6","Site2-WellPN1","Site2-WellPN2"];

m = 1;

for k=1:length(sites)
    site = sites{k}
    
    [T2dist, T2logbins,nmrName] = loadRawNMRdata(site);
    
    % load data file
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    depthsAll = z;
    logT2peakVal = []
    nmrDepths = d(1:end,1);
    
    for j = 1:length(nmrDepths)
        depth = nmrDepths(j)
        T2depths = T2dist(:,1);
        T2ind = (round(T2depths,4) == depth)
        T2interval = T2dist(T2ind,:);
        [T2peakAmp, T2peakInd] = max(T2interval(2:end));
        logT2peakVal(j) = T2logbins(T2peakInd);
    end
    
    %logT2peakVal = logT2peakVal';
    T2peakVal = 10.^(logT2peakVal);
    logT2ML = log10(T2ML);
    logK = log10(Dk);
    
if (site == "Site1-WellG6")


    % Using T2B from Keating and Knight 2007
    % for quartz sand and ferrihydrite-coated sand

    depthCutoff = 5.8;
    T2B = 2.0680;

    Dk = Dk(depthsAll>depthCutoff);
    phi = phi(depthsAll>depthCutoff);
    T2ML = T2ML(depthsAll>depthCutoff);
    SumEch = SumEch(depthsAll>depthCutoff);
    logK = logK(depthsAll>depthCutoff);
    logT2ML = logT2ML(depthsAll>depthCutoff);
    
%     T2B = ones(length(depthsAll),1);
%     T2B(depthsAll > depthCutoff) = 3.048;
%     T2B(depthsAll <= depthCutoff) = 2.433;
elseif (site == "Site1-WellG5")

    % Using T2B from Keating and Knight 2007
    % for quartz sand and ferrihydrite-coated sand

    depthCutoff = 4;
    T2B = 2.0680;

    Dk = Dk(depthsAll>depthCutoff);
    phi = phi(depthsAll>depthCutoff);
    T2ML = T2ML(depthsAll>depthCutoff);
    SumEch = SumEch(depthsAll>depthCutoff);
    logK = logK(depthsAll>depthCutoff);
    logT2ML = logT2ML(depthsAll>depthCutoff);
    
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
     DPP_K{k} = Dk;

    [KGM_lk{k}, bestTau{k}, bestRho{k}, r{k}] = grid_search_kgm(logT2ML, phi, logK, m, T2B);
    
    KGM_K{k} = 10.^KGM_lk{k};
    
    bestError{k} = computeError(DPP_K{k}, KGM_K{k});
    [errorSign{k}, errorFactor{k}] = estimateKdiffFactor_withSign(DPP_K{k},KGM_K{k},1);
    medianErrorFactor(k) = median(errorFactor{k});

end
