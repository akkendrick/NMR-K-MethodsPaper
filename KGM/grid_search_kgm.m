function [KGM_lk, besttau, bestrho, r] = grid_search_kgm(lT2, phi,  kk, m, T2B)
% this function does a grid search over specified parameter values in the
% SDR equation

% put in linear space
T2 = 10.^lT2; 
K = kk; 

% set bounds on A = g/(8*tau^2*rho^2), B = rho^2
% NOTE: TAU and RHO ARE ALREADY IN LOG SPACE
max_log10tauMod = 0; % tau = 1 
min_log10tauMod = -2; % tau = 10

max_log10rho = -2; 
min_log10rho = -6.6; 

% Size of the grid to search
nm = 1000; 

% tau mod is in log space still, rhospace is now in linear space via the
% logspace command 
logTauMod = linspace(min_log10tauMod, max_log10tauMod, nm); 
rhospace = logspace(min_log10rho, max_log10rho, nm); 

parameterMatrixLogTau = logTauMod .* ones(nm,nm);
parameterMatrixRho = rhospace' .* ones(nm,nm);

parameterMatrixLogTau = (fliplr(parameterMatrixLogTau));
parameterMatrixRho = (fliplr(parameterMatrixRho));

% KGM equation
% lK_kgm = @(A,B, T, T2, lphi) (A) + lphi ...
%     + 2*log10( - (D(T)/B) + sqrt((D(T)/B)^2 + ((4*D(T)*Tb(T)*T2)./(Tb(T) - T2)))); 
%lK_kgm = @(A, B, T, T2, phi) KGMfun([A, B], [T2]); 
lK_kgm = @(A, B, C, T, T2, phi, T2B) KGMfun([A, B, C, T], [T2, phi], T2B); 

%Set temperature in the model by estimating average GW temperature at a
%location in Celsius
T = 7;

% run grid search
[r] = deal(zeros(length(logTauMod), length(rhospace)));  
for mloop = 1:length(rhospace)
    for bloop = 1:length(logTauMod)   
        rhoStep = parameterMatrixRho(mloop,bloop);
        logTauModStep = parameterMatrixLogTau(mloop,bloop);
        
        % T2B = NaN means it is estimated using relationship in Dlugosch
        % for T
        lKp = lK_kgm(logTauModStep, rhoStep, m, T, T2, phi, T2B);
        r(mloop, bloop) = norm(lKp - kk); 
    end
end

% % compute best-fitting values
% ind1 = Bspace == 1; 
% ind2 = Bspace == 2; 
% ind4 = Bspace == 4; 
% ind0 = Bspace == 0; 
% 
% r1 = r(ind1, :) == min(r(ind1,:)); 
% r2 = r(ind2, :) == min(r(ind2,:)); 
% r4 = r(ind4,:) == min(r(ind4,:)); 
% r0 = r(ind0,:) == min(r(ind0,:)); 
% 
% bestb(1)= log10(Aspace(r1)); 
% bestb(2)= log10(Aspace(r2)); 
% bestb(3)= log10(Aspace(r4)); 
% bestb(4)= log10(Aspace(r0)); 

s = size(r); 
bestp = min(r(:)); 

indp = find(r(:) == bestp); 
[ii, jj] = ind2sub(s, indp);
% bestA = tauspace(ii); 
% bestB = rhospace(jj); 

bestA = parameterMatrixLogTau(ii,jj);
bestB = parameterMatrixRho(ii,jj);

lkpred = lK_kgm(bestA, bestB, m, T, T2, phi, NaN);

figure; scatter(kk, lkpred), axis equal, hold on, refline(1,0)
xlabel('Observed K')
ylabel('Predicted K')
%title(strcat('KGM K comparison for', site,' m=',string(m)))

taus = 1./sqrt(10.^(logTauMod)); 
tausParamMatrix = 1./sqrt(10.^parameterMatrixLogTau);

% Put rho into log space for plotting
modRhosParamMatrix = log10(parameterMatrixRho);

errorMatrix = log10(r/length(phi));
display('ii')
ii
display('jj')
jj

% Make sure all variables are displayed in LINEAR SPACE
display('tau (ii,jj)')
besttau = tausParamMatrix(ii,jj)
display('rho (ii,jj)')
bestrho = bestB
display('Parameter Error (ii,jj)')
errorMatrix(ii,jj)

KGM_lk = lkpred;
% Plot grid results
figure; 
% subplot(4, 4, [1:3,5:7,9:11, 13:15])
hold on
%imagesc(fliplr(taus), fliplr(log10(rhospace)), errorMatrix)
%imagesc([1 10],[0 -6.6],flipud(fliplr(errorMatrix)))
%imagesc([1 10],[0 -6.6],errorMatrix)

%uimagesc(fliplr(taus),fliplr(log10(rhospace)),errorMatrix)
uimagesc(fliplr(taus),(log10(rhospace)),errorMatrix)
box on
grid on

% tausParamMatrix = fliplr(tausParamMatrix);
% rhosParamMatrix = flipud(rhosParamMatrix);
% errorMatrix = flipud(fliplr(errorMatrix));

bestp = min(errorMatrix(:)); 

indp = find(errorMatrix(:) == bestp); 
[ii, jj] = ind2sub(s, indp);

plot(tausParamMatrix(ii,jj), modRhosParamMatrix(ii,jj), '*w')

%% Plot grid search results

caxis([-1.2,0])
% plot( bestb, [1, 2, 4, 0], 'sw')
colorbar('Ticks', [-1.2:.2:0])
% colorbar
colormap(jet)
xlim([min((taus)), max((taus))])
ylim([min(log10(rhospace)), max(log10(rhospace))])
ylabel('log_{10} (\rho)')
xlabel('\tau')
set(gca,'FontSize',14)

%str = [strcat(site, ' m =',string(m)),...
%    strcat('\tau =', string(besttau), ' \rho =', string(bestrho))];
%title(str)

% subplot(4,4,16)
% hold on
% plot(log10(Aspace), r(ind0,:))
% title('(e)  m = 0')
% xlabel('log_{10}(b)')
% ylabel('Residual')
% 
% subplot(4,4,12)
% hold on
% plot(log10(Aspace), r(ind1,:))
% title('(d)  m = 1')
% xlabel('log_{10}(b)')
% ylabel('Residual')
% 
% subplot(4,4,8)
% hold on
% plot(log10(Aspace), r(ind2,:))
% title('(c)  m = 2')
% xlabel('log_{10}(b)')
% ylabel('Residual')
% 
% subplot(4,4,4)
% hold on
% plot(log10(Aspace), r(ind4,:))
% title('(b)  m = 4')
% xlabel('log_{10}(b)')
% ylabel('Residual')
% 
% set(gcf,'OuterPosition',[100 100 1000 800], 'PaperPositionMode','auto','color',[1 1 1])

% % output best b values
% bs = [Aspace(ii); bestb(:)]; 
% ms = [Bspace(jj); [0,1, 2, 4]']; 
end