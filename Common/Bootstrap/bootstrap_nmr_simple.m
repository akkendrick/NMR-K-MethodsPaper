function [b] = bootstrap_nmr_simple (data, n, m)
%Function called by bootstrap_fixed_m_n to run the bootstrap
% Data is [logT2ML, logPhi, logK]

% y is log K
y = data(:,end);

% x is [ones, logT2ML, logPhi]
% ones is for solving for b, SDR model is b*phi^m*T2ML^n, this is
% similar to these types of equations just in log form to weight K
% variations equally for the LLS model
x = [ones(size(data,1),1), data(:,1), data(:,2)]; 

ft = fittype('poly11');
cterm1 = m*x(:,3);
cterm2 = n*x(:,2);
x = [cterm1, cterm2]; 
fo = fitoptions('method','LinearLeastSquares','Lower',[-Inf, 1, 1],'Upper',[Inf, 1, 1]);
model = fit(x,y,ft,fo);
b(1) = model.p00;
b(2) = n;
b(3) = m;
       
end
