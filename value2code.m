function [yr, yc, Xr, Xc] = value2code(N, Mr, Mc, y, phi, sd_phi)
%%       Convert Component Values to Parameter Matrix X 

% K :   the number of used components
% Mr :  the number of standard components of resistor
% Mc :  the number of standard components of capacitor
% w :   vector of frequency sample points
% y :   vector of component values
% phi: the vector of standard components
% sd_phi: the variance of standard components

% Copyright (c) 2018, Kunpeng Wang.
% Email: wkphnzk@163.com
%

Xr = zeros(N,Mr);
Xc = zeros(N,Mc);
for i=1:N
    for j=1:Mr
        r_err = max(y(i) - Xr(i,1:j)*phi(1:j),0);
        Xr(i,j) = floor(r_err/phi(j));
        uncertain = sum(Xr(i,1:j)*sd_phi(1:j)); 
        r_err = max(y(i) - Xr(i,1:j)*phi(1:j),0);
        if uncertain > r_err    % if the uncertainty > error stop
            break;
        end
    end

    for j=1:Mc
        c_err = max(y(N+i) - Xc(i,1:j)*phi(Mr+1:Mr+j),0);
        Xc(i,j) = floor(c_err/phi(Mr+j));
        uncertain = sum(Xc(i,1:j)*sd_phi(Mr+1:Mr+j)); 
        c_err = max(y(N+i) - Xc(i,1:j)*phi(Mr+1:Mr+j),0);
        if uncertain > c_err    % if the uncertainty > error stop
            break;
        end
    end
end

yr = Xr*phi(1:Mr);
yc = Xc*phi(Mr+1:end);
