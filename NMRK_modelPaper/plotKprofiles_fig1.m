% Plot all K, T2ml, water table profiles

clear
%close all

m = 1;
n = 2;


figureson = 0;
wDirect = 0;
saveData = 0;
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

%load bProfileData.mat
%load bestKModels.mat

sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};
figNames = {'Adams G5','Adams G6','Plainfield Lake PN1','Plainfield Lake PN2'}
 
gammaNames = {'Site1-G5-Well1-gamma-EMI-bLS.csv','Site1-G6-Well2-gamma-EMI-bLS.csv','Site2-Well1-gamma-EMI-bLS.csv','Site2-Well2-gamma-EMI-bLS.csv'};

%waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface %NOTE G5/G6 cased below clay,
%   need to use nearby water level from above the clay, for G5 + G6 using
%   water level from well well G2 cased above the New Rome Clay (rel to gs)
waterTable = [2.0469,2.1248,5.0285,4.7476];
plotNum = 1;
f1 = figure('Renderer', 'painters', 'Position', [10 10 1000 1150])

 for jj = 1:length(sites)
    baseName = sites{jj}
    figName = figNames{jj}
    % Add SDR K estimates
    
    
    
%      baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
%      gammaBaseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/WI_gamma-EMI-bLS_csvfiles';

    baseDir = 'C:\Users\kenta\Dropbox\Research\Alex-Rosemary\Papers\Kmodel_comparision_paper\howTo\Common\Field Data\USGS Data\';
    gammaBaseDir = 'C:\Users\kenta\Dropbox\Research\Alex-Rosemary\Papers\Kmodel_comparision_paper\howTo\Common\Field Data\USGS Data\WI_gamma-EMI-bLS_csvfiles';
    
    [T2dist,T2logbins,SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
        oneDVectorsUniform, nmrName] = loadAllRawNMRdata(baseName);

    in3 = [baseDir baseName '/' strcat(baseName,'_DPP_filt.txt')];
    in7 = [gammaBaseDir '/' gammaNames{jj}];
    
    DPPdat = load(in3); 
    gammaEMIdata = csvread(in7,3,0);
    
%     T2depths = T2dist(:,1);
%     T2dist = T2dist(:,2:end);
    [K,z,k_btsrp,k_mcmc,k_direct,bestFitMatrix,b_boot,totalErrorEstimate] = computeProfile_modT2ML(baseName,n,m,figureson,wDirect,saveData);
 
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(baseName); 
    
    for j = 1:length(b_boot)
        SDR_K_all{j} = SDR_K(b_boot(j),m,n,phi,T2ML);
    end
    
    SDR_K_best = SDR_K(median(b_boot), m,n,phi,T2ML);

    NMRphi = oneDVectors(:,2);

    Dk = DPPdat(:,2)*1.16e-5; % converts K from m/day to m/s
    z_dk = DPPdat(:,1);

    origT2dist = T2dist;
    
    %%
    zNMR = z';
   
    %%
    
    gammaEMIoffset = [0,0,0.75, 0.75];
    
    if jj == 1 || jj == 2
        coreSampleLoc = [9 14 17.5 22.5 26 27.5 31 36 39 42.5 47 51.73 54 59 63]*0.3048;
        print jj
    else
        coreSampleLoc = [];
    end
        
    plotWellProfiles_NMR_DPP(waterTable(jj),z,T2dist,T2logbins,DPPdat,SDR_K_all, SDR_K_best, plotNum, f1, coreSampleLoc)
    title(figName)
    
    
%     fileName = strcat(baseName,'_profile.fig');
%     savefig(fileName)
%     
%     fileName = strcat(baseName,'_profile.png');
%     print('-dpng','-r300',fileName)
%     
    
    %title(baseName)
    plotNum = plotNum + 2;

    
 end
 
fileName = 'allSites_Kline_profile.png';
print('-dpng','-r300',fileName)

    