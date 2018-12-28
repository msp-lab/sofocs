function [UC, D_sigma] = uncertainty_eval(N, structure, kappa, q, ...
                                          delta, Xr, Xc, ItNum, stdrc_dir)
%% Monte Carlo based Unvertainty Estimation Algorithm
%% Input:
%  N: system order of franctance 
%  structure : The structure of fractance, Desired/Chain/Tree/Ladder/ZPK
%  mag_true : vector of desired Magnitude Response Curve
%  q : the order of fractional-order chaos system
%  delta : dB, the maximum discrepancy between desired H(s) and \hat{H}(s)
%  Xr: parameter matrix of resistor
%  Xc: parameter matrix of capicitor
%  ItNum: repeat times
%% Output:
%  UC: uncertainty in circuit implementation
%  D_sigma: maximum discrepancy between \hat{H}(s) and H(s)

%
% Copyright (c) 2018, Kunpeng Wang.
% Email: wkphnzk@163.com
%

stdrc = load('stdrc.mat', ...
             'rstd', 'cstd', 'sd_rstd', 'sd_cstd', 'omega');
omega = stdrc.omega;
rstd = stdrc.rstd;
cstd = stdrc.cstd;
sd_rstd = stdrc.sd_rstd;
sd_cstd = stdrc.sd_cstd;

mag_true = mag_curve('Desired', q, delta, omega);

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

disp(['++-- Fractional order q=' num2str(q)])
% rstd(j), cstd(j): Mean of the non-truncated Gaussian variables
% sd_rstd(j), sd_cstd(j): Standard deviation of the non-truncated Gaussian
disp('Start uncertainty analysis:')
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
    mag_opt = mag_curve(structure, q, delta, omega, [yr; yc; kappa], N);
    AvdB = 20*log10(mag_true);
    nAvdB = 20*log10(mag_opt);
    D_sigma(i) = max(abs(nAvdB - AvdB));
end
     
UC = sum(D_sigma>delta) / ItNum;

disp(['E{D_sigma}=' num2str(mean(D_sigma)) ', '...
      'delta=' num2str(delta) ', ' ...
      'Uncertainty UC=' num2str(UC)]);





end