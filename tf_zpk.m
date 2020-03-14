function [Hs, tf_N] = tf_zpk(q, delta, omega)
%%       Transfer Function H(s) of the "pole/zero" Method

%% Input
% q : the order of fractional-order chaos system
% delta : dB, the maximum discrepancy between desired H(s) and \hat{H}(s)
% omega :  vector of frequency sample points

%% Output
% Hs: the transfer function of the "pole/zero" Method
% tf_N: the system order of H(s)

% Copyright (c) 2020, Kunpeng Wang.
% 
    tau0 = 100;  % the relaxation time consistant
    p_0 = 1/tau0*10^(delta/(20*q));
    a = 10^(delta/(10*(1-q)));
    b = 10^(delta/(10*q));
    w_max = max(omega)/(2*pi);
    % the total number of the poles in H(s)
    tf_N = 1 + floor(log10(w_max/p_0)/log10(a*b)); 

    z = zeros(1, tf_N);
    p = zeros(1, tf_N+1);
    for i=0:tf_N
       if i < tf_N
            z(i+1) = ((a*b)^i*a*p_0); % zeros
       end
       p(i+1) = ((a*b)^i*p_0);        % poles
    end
    gain = b^tf_N * p_0;  % gain

    Hs = zpk(z, p, gain); % transfer function H(s)
    tf_N = tf_N + 1;      % the order of transfer function
end

