% Calculate K from known rho and tau values
clear

siteList = ["Site1-WellG5","Site1-WellG6","Site2-WellPN1","Site2-WellPN2"];

n = 2;
m = 1;

lK_kgm = @(A, B, C, T, T2, phi, T2B) KGMfun([A, B, C, T], [T2, phi], T2B); 
T = 7;

%rho = 8.52E-05; % Value from Wyoming Ren et al. 2018
%rho = 19.01E-06; % Max Wisconsin core value 
rho = 0.0003876 ;% Best Calibrated value

%tau = 1.38; % Value from Wyoming Ren et al. 2018
tau = 1.7783; % Best Calibrated value

for i = 1:length(siteList)
    siteName = siteList{i};
    [T2dist, T2logbins, nmrName{i}] = loadRawNMRdata(siteName);
    
    [d, DPP_K{i}, T2ML, phi{i}, z{i}, SumEch{i},K_SOE{i}, logK{i}, logT2ML{i}, logPhi{i}, SumEch_3s{i}, SumEch_twm{i}, ...
    SumEch_twm_3s{i}] = loadnmrdata2_Ksoe(nmrName{i}); 

    % Compute the maximum T2

    logT2peakVal = []
    nmrDepths = d(:,1);
    for j = 1:length(nmrDepths)
        depth = nmrDepths(j)
        T2depths = T2dist(:,1);
        T2ind = (round(T2depths,4) == depth)
        T2interval = T2dist(T2ind,:);
        [T2peakAmp, T2peakInd] = max(T2interval(2:end));
        logT2peakVal(j) = T2logbins(T2peakInd);
    end
    
    logT2peakVal_step{i} = logT2peakVal';
    T2peakVal{i} = 10.^(logT2peakVal');
    logK{i} = log10(DPP_K{i});
    
    if (siteName == "Site1-WellG6")


        % Using T2B from Keating and Knight 2007
        % for quartz sand and ferrihydrite-coated sand

        depthCutoff = 5.8;
        T2B = ones(length(nmrDepths),1);
        T2B(nmrDepths > depthCutoff) = 3.048;
        T2B(nmrDepths <= depthCutoff) = 2.433;

    elseif (siteName == "Site1-WellG5")

        % Using T2B from Keating and Knight 2007
        % for quartz sand and ferrihydrite-coated sand

        depthCutoff = 4;
        T2B = ones(length(nmrDepths),1);
        T2B(nmrDepths > depthCutoff) = 3.048;
        T2B(nmrDepths <= depthCutoff) = 2.433;
    else
        T2B = ones(length(nmrDepths),1);

        T2B = T2B .* 3.048;

        %T2B = 2.0680; % Fix T2B for other sites where we don't have data
        % this is T2B from KGM model equation for T2B at 7 deg C estimated from
        % map of groundwater temperature
        % (https://pubs.usgs.gov/wsp/0520f/report.pdf) 
        % "Temperature of water available for industrial use in the United States: 
        % Chapter F in Contributions to the hydrology of the United States,
        % 1923-1924" by W.D. Collins
        % Equation for T2B from Dlugosch et al. 2013
    end
    T2B_site{i} = T2B;
end

T2B_all = vertcat(T2B_site{:});
logT2peakVal_all = vertcat(logT2peakVal_step{:});
T2peakVal_all = 10.^logT2peakVal_all;
phi_all = vertcat(phi{:});
logK_all = vertcat(logK{:});
DPP_K_all = vertcat(DPP_K{:});

rhoVals = ones(size(T2peakVal)).*rho;
tauVals = ones(size(T2peakVal)).*tau;

    
lKp1 = lK_kgm(log10(tau), rho, m, T, T2peakVal_all, phi_all, T2B_all);    
KGM_K = 10.^lKp1;

bestError = computeError(DPP_K_all, KGM_K);
[errorSign, errorFactor] = estimateKdiffFactor_withSign(DPP_K_all,KGM_K,1);
medianErrorFactor = median(errorFactor)



%% Plot results
totalKGMK = KGM_K;
totalDPPK = DPP_K_all;

figure(1)

hold on
    
grid on 
box on

scatter(totalDPPK, totalKGMK,60,'Filled')

plot(totalDPPK,totalDPPK,'k','LineWidth',2,'HandleVisibility','off')
plot(totalDPPK,totalDPPK*10,'k--','HandleVisibility','off')
plot(totalDPPK,totalDPPK*0.1,'k--','HandleVisibility','off')

legend('KGM','Location','northwest')
ylim([5*10^-8,1*10^0])
xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')