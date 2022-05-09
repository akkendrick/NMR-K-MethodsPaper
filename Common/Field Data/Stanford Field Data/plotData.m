% Plot well data together to identify relevant intervals 

load AdamsDPP.mat
load sonicCoreIronData.mat
load sonicCoreT2B.mat

%goodDepths = [2.743 4.267, 5.334, 6.401, 7.925, 8.382, 9.449,...
%    10.97, 11.89, 13.41, 15.77, 16.46, 17.98, 19.51]

%goodDepthsFeet = [9, 14, 17.5, 21, 26, 27.5, 31, 36, 39, 44, 51.73, 54, 59, 64];

newDepthsFeet = [9, 14, 17.5, 22.5, 26, 27.5, 31, 36, 39, 42.5,...
    47, 51.73, 54, 59, 63];

% change 21 ft to one below it
% change 44 ft to 42.5 ft (correct to m) 

% think about 46 ft could take none in this interval so far

% think about how useful 59 ft is
% change 64 to 63 ft

% Due to concern about pref flow pathways, will take 4 inches of core above
% the locations sampled for water

%goodDepths = goodDepths .* 3.28;
goodDepths = newDepthsFeet;
sampleWidth = 0.1016 * 3.28;

goodDepthStart = goodDepths;
goodDepthEnd = goodDepths + sampleWidth;

lineX = [zeros(1,length(goodDepths)); ones(1,length(goodDepths)).*3000];

lineY = [goodDepths; goodDepths];
lineY_extent = [goodDepths + sampleWidth; goodDepths+sampleWidth]; 

figure(1)
subplot(1,3,1)
hold on
plot(AdamsDPP.SouthK,AdamsDPP.SouthDepth.*3.28,'LineWidth',2)
plot(AdamsDPP.NorthK,AdamsDPP.NorthDepth.*3.28,'LineWidth',2)
xlabel('K (m/day)')
title('KGS DPP Data')
ylabel('Depth (ft)')
set(gca,'Ydir','reverse')
xlim([0,60])
ylim([0,70])
line(lineX,lineY,'Color','k')
line(lineX,lineY_extent,'Color','r')

subplot(1,3,2)
hold on
plot(sonicCoreIronData.TotalFe,sonicCoreIronData.Depthm.*3.28,'*')
plot(sonicCoreIronData.Fe2,sonicCoreIronData.Depthm.*3.28,'*')
line(lineX,lineY,'Color','k')
line(lineX,lineY_extent,'Color','r')
title('UV-Vis Data')
xlabel('Concentration uM')

set(gca,'Ydir','reverse')
xlim([0,90])
ylim([0,70])

subplot(1,3,3)
hold on
plot(sonicCoreT2BData.T2Bpeak,sonicCoreT2BData.Depthm.*3.28,'*')
set(gca,'Ydir','reverse')
xlim([0,2500])
ylim([0,70])
line(lineX,lineY,'Color','k')
line(lineX,lineY_extent,'Color','r')
title('NMR T_{2B} Data')
xlabel('T_{2B} (ms)')

