function [logKpred] = KGMfun(x, xdata, T2B)
% this function computes the likelihood for Kpredicted given K_measure

    [aa,bb] = size(x); 
    if aa ~=1 & bb~=1
        if bb > aa
            x = x'; 
        end
    end
    T2 = xdata(:,1);
    if size(xdata,2)==2
        phi = xdata(:,2);        
        lphi = log10(phi(:)); 
    end

    % Hardcoding x 10/5 to include addition of Temp
    Temp = x(:,4);
    m = x(:,3);
    logTauMod = x(:,1);
    rhoNMR = x(:,2);
    
    rho = @(Tt) 1000*(1 - ((Tt+288.94)./(508929*(Tt+68.12))).*(Tt-3.98).^2); % kg/m^3
    eta = @(Tt) 0.0013 - 1.7e-5*Tt;         % Pa -s
    %Tb = @(Tt) 2.2293;
    D = @(Tt) (1.0413 + 0.039828*Tt + 0.00040318*Tt.^2).*1e-9;  % m^2/s 
    g = 9.8;    %m/s^2
    %tort = 1/(1.5^2); 
    
    %Temp = 7.4;  % temperature in degress C 

    if isnan(T2B)
        Tb = @(Tt) 3.3 + 0.044*(Tt - 35);  % seconds
        
        t1 = (rho(Temp)*g)/(8*eta(Temp)); % 

        num2 = 4*D(Temp)*Tb(Temp)*T2;
        denom2 = Tb(Temp) - T2; 

        f12 = (D(Temp)./rhoNMR);  
        SQterm = sqrt(f12.^2 + (num2./denom2)); 
    else
        Tb = T2B;
        
        t1 = (rho(Temp)*g)/(8*eta(Temp)); % 

        num2 = 4*D(Temp).*Tb.*T2;
        denom2 = Tb - T2; 

        f12 = (D(Temp)./rhoNMR);  
        SQterm = sqrt(f12.^2 + (num2./denom2)); 
        
    end
    



%     if size(x,2) ==3
%         ltinv = x(:,1); 
%         B = x(:,2); 
%         m = x(:,3);
%     elseif size(x,2) == 2
%         ltinv = x(:,1); 
%         B = x(:,2); 
%         m = 0; 
%     else
%         ltinv = tort;
%         B = x(1); 
%         m = 0; 
%     end



    
 

    % predicted data
    % Note logTauMod is log10(1/tau^2)
    if size(xdata,2)==2 % include porosity term
        logKpred = logTauMod + log10(t1) + m*lphi + 2*log10(SQterm-f12);
    else %m = 0 case, don't include porosity
        logKpred = ltinv + log10(t1) + 2*log10(SQterm-f12);
    end

end