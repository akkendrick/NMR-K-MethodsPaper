% Run Seevers Models

% Range over pairs of m and n values
close all
clear

siteList = [{'Site1-WellG5'} {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
%siteList = [{'Site1-WellG6'}];

% 
%  siteList = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%    'dpnmr_leque_east','dpnmr_leque_west'};

%  siteList = {'dpnmr_larned_east','dpnmr_larned_lwph',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%    'dpnmr_leque_east','dpnmr_leque_west'};

% %m = [0 1 2 4 0 1 2 4];
%n = [1 1 1 1 2 2 2 2];

m = [1];
n = [2];

%m = [1]
%n = [1]


figureson = 1;
wDirect = 0;
currentRow = -3;

mcmcMatrix = [];
bootMatrix = [];
directMatrix = [];

% % for k = 1:length(siteList)
% % 
% %     computeSeevers(siteList{k},n,m,figureson,0)
% %     title(siteList{k})
% % 
% % end
totalbMatrix = zeros(3,length(m),length(siteList));
totalErrorMatrix = zeros(3,length(m),length(siteList));

tic

if isempty(m) && isempty(n)
    currentFitMatrix = [];
    currentErrorMatrix = [];
    for i = 1:length(siteList)
        siteName = siteList{i}

        [K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,b_boot,...
            Seevers_bs{j},totalErrorEstimate] = computeSeevers(siteName,[],[],figureson,wDirect);

        currentFitMatrix = [currentFitMatrix bestFitMatrix];
        currentErrorMatrix = [currentErrorMatrix totalErrorEstimate];
        
        b_boot_all{i} = b_boot;

    end
else
    for j = 1:length(m)
        currentRow = currentRow + 4;
        currentFitMatrix = [];
        currentErrorMatrix = [];
        
        tempb = zeros(3,length(siteList));
        tempError = zeros(3,length(siteList));

        for i = 1:length(siteList)
            siteName = siteList{i}

            [K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,...
                b_boot,Seevers_bs{j},totalErrorEstimate] = computeSeevers(siteName,n(j),m(j),figureson,wDirect);

            tempb(:,i) = bestFitMatrix(1,:)';
            tempError(:,i) = totalErrorEstimate(1,:)';
            currentFitMatrix = [currentFitMatrix bestFitMatrix];
            
             b_boot_all{i} = b_boot;


        end

        matrixKey(1,j) = n(j);
        matrixKey(2,j) = m(j);
        
        totalbMatrix(:,j,:) = tempb;
        totalErrorMatrix(:,j,:) = tempError;
    end
end

toc

save('Seevers_dat.mat','siteList','m','n','matrixKey','totalbMatrix','b_boot_all','totalErrorMatrix','currentFitMatrix')

%%
