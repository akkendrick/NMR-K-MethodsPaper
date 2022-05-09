% Try calculating SOE from raw data
close all
clear


site = {'Site1-WellG5'}
figureson = 1; 

[T2dist,T2logbins,SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
    oneDVectorsUniform, nmrName] = loadAllRawNMRdata(site);

depths = T2dist(:,1);
SEdecay = SEdecay(:,2:end);
SEdecayUniform = SEdecayUniform(:,2:end);
noise = oneDVectors(:,14);
porosity = oneDVectors(:,2);
T2dist = T2dist(:,2:end);

T2linbins = 10.^T2logbins;

%% 
slice = 30;
smoothing = 50;

figure(1)
hold on
plot(SEdecayTime, SEdecay(slice,:))
%%
testSOE = sum(SEdecay(slice,:))