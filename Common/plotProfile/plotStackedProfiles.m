% Plot basic K, T2ml, water table, phi profiles

clear
%close all

%load bProfileData.mat
%load bestKModels.mat
sites = {'Site1-WellG5','Site1-WellG6'};
 
 
gammaNames = {'Site1-G5-Well1-gamma-EMI-bLS.csv','Site1-G6-Well2-gamma-EMI-bLS.csv','Site2-Well1-gamma-EMI-bLS.csv','Site2-Well2-gamma-EMI-bLS.csv'};

%waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface %NOTE G5/G6 cased below clay,
%   need to use nearby water level from above the clay, for G5 + G6 using
%   water level from well well G2 cased above the New Rome Clay (rel to gs)
waterTable = [2.0469,2.1248,5.0285,4.7476];
    
 for jj = 1:length(sites)
    baseName = sites{jj}
%      baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
%      gammaBaseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/WI_gamma-EMI-bLS_csvfiles';

    baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';
    gammaBaseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\WI_gamma-EMI-bLS_csvfiles';
    
    [T2dist{jj},T2logbins{jj},SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
        oneDVectorsUniform, nmrName] = loadAllRawNMRdata(baseName);

    in3 = [baseDir baseName '/' strcat(baseName,'_DPP_filt.txt')];
    in7 = [gammaBaseDir '/' gammaNames{jj}];
    
    DPPdat{jj} = load(in3); 
    gammaEMIdata{jj} = csvread(in7,3,0);
    
%     T2depths = T2dist(:,1);
%     T2dist = T2dist(:,2:end);

    [d, K, T2ML{jj}, phi, z{jj}, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(baseName); 

    NMRphi{jj} = oneDVectors(:,2);

    Dk{jj} = DPPdat{jj}(:,2)*1.16e-5; % converts K from m/day to m/s
    z_dk{jj} = DPPdat{jj}(:,1);

 end
   
    %%
    
    gammaEMIoffset = [0,0,0.75, 0.75];
     
%         plotWellProfiles_wLith(K,NMRphi,waterTable(kk),z,T2ML,T2dist,T2logbins,[],{'K DPP'},{'*'},bProfile{kk},...
%             lithStart{kk}, lithEnd{kk},lithUnits{kk},dlubacModel{kk},bestSDRModel{kk},meanT2B)
    
    figure('Renderer', 'painters', 'Position', [10 10 300 800])

    plotWellProfiles_stacked(sites,DPPdat)


%     fileName = strcat(baseName,'_profile.fig');
%     savefig(fileName)
%     
%     fileName = strcat(baseName,'_profile.png');
%     print('-dpng','-r300',fileName)
%     
%     fileName = strcat(baseName,'_profile.svg');
%     print('-dsvg','-r300',fileName)
    
    %title(baseName)
    
 
    