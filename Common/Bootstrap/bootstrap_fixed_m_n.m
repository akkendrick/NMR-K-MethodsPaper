function [bs, ns, ms, logb] = bootstrap_fixed_m_n (data, Nboot, n, m)

% call function, depending on whether or not m and n are specified
bhat = bootstrp(Nboot, @(x,y,z) bootstrap_nmr_simple(x,n,m), data); % m = 0, n = 2

bs = bhat(:,1); 
ns = bhat(:,2);
ms = bhat(:,3); 


% transform back into linear space
logb = bs; 
bs = 10.^logb;

% m and n should be the same

end