% Range over pairs of m and n values
close all
clear

%Wisconsin Data
siteList = [{'Site1-WellG5'}...
   {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
    
% Maurer and Knight
%  siteList = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%    'dpnmr_leque_east','dpnmr_leque_west'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define fitting variables and K range
m = 0;
n = 2;
Kcutoff = 5*10^-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureson = 0;
wDirect = 0;
currentRow = -3;
saveData = 0;

        
for i = 1:length(siteList)
    siteName = siteList{i};
    disp(strcat('On site: ', siteName))
    tic
    
       [K,z,b_boot,k_boot,totalErrorEstimate] = computeProfile_Kfilt(siteName,...
            n,m,Kcutoff,figureson,wDirect,saveData);

    b_median{i} = median(b_boot);
    b_all{i} = b_boot;
    m_site{i} = m;
    n_site{i} = n;
    NMR_K{i} = k_boot;
    DPP_K{:,i} = K;
    
    if isnan(NMR_K{i})
        allErrorFactors{i} = NaN;
        medianErrorFactor{i} = NaN;
    else
        allErrorFactors{i} = estimateKdiffFactor(DPP_K{i},NMR_K{i},1);
        medianErrorFactor{i} = median(allErrorFactors{i});
    end
    
    toc
        

end

        
%save('SDR_T2lm_Kfilt_bestFit_092821_m1_n2.mat','siteList','m','n','matrixKey','totalbMatrix','totalmMatrix','totalErrorMatrix','totalKMatrix')
            
            