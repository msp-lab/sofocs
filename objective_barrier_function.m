function G = objective_barrier_function(ystate, stage, q, delta,...
                                omega, structure, mag_true, phase_true,...
                                N, Mr, Mc, phi, sd_phi)
%%       Objective Function for the Fractional-order Chaos 
%%              Circuit Parameter Selection Problem

% q : the order of fractional-order chaos system
% delta : dB, the maximum discrepancy between desired H(s) and \hat{H}(s)
% omega : vector of frequency sample points
% structure : The structure of fractance, Desired/Chain/Tree/Ladder/ZPK
% mag_true : vector of desired Magnitude Response Curve
% phase_true: vector of desired phase of Magnitude Response Curve
% delta : dB, the maximum discrepancy between desired H(s) and \hat{H}(s)
% K :  the number of used components
% Mr :  the number of standard components of resistor
% Mc :  the number of standard components of capacitor
% phi: the vector of standard components
% sd_phi: the variance of standard components

% Copyright (c) 2018, Kunpeng Wang.
% Email: wkphnzk@163.com
%


assert(ischar(stage)&&isrow(stage),'stage must be a 1xN char.')

kappa = ystate(end);
y = ystate(1:end-1)';

lambda = 0.5;
lambda1 = 100;  % trade-off between sparseness and uncertainty
lambda2 = 0.65; %trade-off between ampitude and phase

switch stage
    case 'INIT'
        % Compute new response curve 
        [mag_new, phase_new] = mag_curve(structure, q, delta, omega,...
                                         ystate, N);
        AvdB = 20*log10(mag_true);
        nAvdB = 20*log10(mag_new);

        % Compute the loss function of response curve error
        err_curve = (1-lambda2)*sum(abs(nAvdB - AvdB)) + ...
                    lambda2*sum(abs(phase_new - phase_true));
        G = err_curve;
        
    case 'ITER'
        % Compute component values sparse representation code
        [yr, yc, Xr, Xc] = value2code(N,Mr,Mc,y,phi,sd_phi);
        num_stdcomp = lambda * sum(Xr(:)) / N +...
                      (1 - lambda) * sum(Xc(:)) / N;
        % Compute the ralitive uncertainty of y
%         sd_r = Xr * sd_phi(1:Mr);
%         sd_c = Xc * sd_phi(Mr+1:end);
%         sd = [sd_r; sd_c];  % sigma_phi uncertainty
%         yv = [yr; yc];
%         cv = sd ./ yv;
%         uncertain = sum(cv(:));
        ItNum = 300;
        UC = uncertainty(N, structure, delta, q, omega, Mr, Mc, Xr, Xc, ...
                         kappa, mag_true, phi, sd_phi, ItNum);
        uncertain = UC;
        
        ystate_new = [yr', yc', kappa];
        % Compute new response curve 
        mag_new = mag_curve(structure, q, delta, omega, ystate_new, N);
        AvdB = 20*log10(mag_true);
        nAvdB = 20*log10(mag_new);
        dist = abs(nAvdB - AvdB);
        % Calculate nonlinear contraint
        max_d = max(dist);
        p = 1.01;
        if max_d > delta
            noncon = log(p*max_d/delta)/log(p) + 1;
        else
            noncon = (max_d/delta).^2; 
        end
        G =  num_stdcomp + lambda1 * uncertain + noncon;
        
%         disp(['G=', num2str(G)])
%         disp(['  max_d=', num2str(max_d)])
%         disp(['  noncon=', num2str(noncon)])
%         disp(['  num_stdcomp=' num2str(num_stdcomp)])
%         disp(['  uncertain=' num2str(lambda1 * uncertain)])

    otherwise 
        error('Stage "%s" is not supported.',stage)
end

end

% Sub-function of uncertainty
function UC = uncertainty(N, structure, delta, q, omega, Mr, Mc, Xr, Xc,...
                          kappa, mag_true, phi, sd_phi, ItNum)
    rstd = phi(1:Mr);
    cstd = phi(Mr+1:Mr+Mc);
    sd_rstd = sd_phi(1:Mr);
    sd_cstd = sd_phi(Mr+1:Mr+Mc);
    ar = rstd - 3*sd_rstd;  % Left bound of resistor
    br = rstd + 3*sd_rstd;  % Right bound of resistor
    ac = cstd - 3*sd_cstd;  % Left bound of capacitor
    bc = cstd + 3*sd_cstd;  % Right bound of capacitor

    %% Generate a big number of variables
    % Number of RV to generate: ItNum
    % Random variable generation
    rstd_v = rstd;
    cstd_v = cstd;
    if min(size(Xr))>1
        r_mask =  sum(Xr);
    else
        r_mask =  Xr;
    end

    if min(size(Xc))>1
        c_mask =  sum(Xc);
    else
        c_mask =  Xc;
    end

%     disp(['++-- Fractional order q=' num2str(q)])
    % rstd(j), cstd(j): Mean of the non-truncated Gaussian variables
    % sd_rstd(j), sd_cstd(j): Standard deviation of the non-truncated Gaussian
%     disp('Start uncertainty analysis:')
    D_sigma = zeros(1, ItNum);
    for i = 1:ItNum
        for j = 1:length(rstd)
            if r_mask(j) > 0
                rstd_v(j) = rtnorm(ar(j),br(j),rstd(j),sd_rstd(j));
            end
        end

        for j = 1:length(cstd)
            if c_mask(j) > 0
                cstd_v(j) = rtnorm(ac(j),bc(j),cstd(j),sd_cstd(j));
            end
        end

        % Mean of sum of the non-truncated Gaussian variables
        yr = Xr * rstd_v;
        yc = Xc * cstd_v;

        % get the system order of transfer function
        mag_opt = mag_curve(structure, q, delta, omega, ...
                            [yr; yc; kappa], N);
        AvdB = 20*log10(mag_true);
        nAvdB = 20*log10(mag_opt);
        D_sigma(i) = max(abs(nAvdB - AvdB));
    end

    UC = sum(D_sigma>delta) / ItNum;

end

