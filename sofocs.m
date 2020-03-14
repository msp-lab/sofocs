function [Xr_zpk, Xc_zpk, yr_zpk, yc_zpk, kappa_zpk, ...
          kappa, Xr, Xc, yr, yc, D_sigma] = sofocs(N, omega, struture,...
                                       mag_true, phase_true,...
                                       mag_zpk, phase_zpk,...
                                       Mr, Mc, q, delta, phi, sd_phi, ...
                                       option)
%% Finding an Approximation Optimal Solution for the Franctance 
%% in Chaotic Circuit System

%% Input:
%  N: system order of franctance 
%  omega: vector of frequency sample points
%  struture : The structure of fractance, Desired/Chain/Tree/Ladder/ZPK
%  mag_true : vector of desired Magnitude Response Curve
%  phase_true: vector of desired phase of Magnitude Response Curve
%  mag_zpk :  Magnitude Response Curve of zpk method
%  phase_zpk: phase of Magnitude Response Curve of zpk method
%  q : the order of fractional-order chaos system
%  delta : dB, the maximum discrepancy between desired H(s) and \hat{H}(s)
%  phi: the vector of standard components
%  sd_phi: the variance of standard components
%  option: 'PoleZero', 'Proposed'
%% Output:
%  X: parameter matrix of our method
%  X_zpk: parameter matrix of "pole/zero" method
%  y: component values of our method
%  y_zpk: component values of "pole/zero" method
%  kappa: system gain adjustment factor
%  gain: system gain G
%  Hs_rc: RC transfer function \hat{H}(s)
%  D_sigma: maximum discrepancy between \hat{H}(s) and H(s)

%
%  Copyright (c) 2020, Kunpeng Wang.
% 

%% Start Sparse Optimization
disp('---- Start Sparse Optimization')
%% Set Initial Values of Resistors and Capicitors in Fractance
disp('---->> Set Initial Values of Components')
% For example, we use the combined RC fractance structure
K = 2*N; % the number of components

% Bounds
lb_x = min(phi) * ones(K,1);
lb_k = 1e-8; 
lb = [lb_x; lb_k];

ub_x =  max(phi)* ones(K,1);
ub_k = 1e-5;
ub = [ub_x; ub_k];
init_yr = 1e5*ones(N,1);
init_yc = 1e4*ones(N,1);
init_kappa = 2e-7;
init_ystate = [init_yr', init_yc', init_kappa];

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%          Optimization for Finding Intial Point
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Create Handle to Custom Output Function
init_output = @(a1,a2,a3)optim_plot(a1, a2, a3, omega, struture,...
                          mag_true, phase_true, mag_zpk, phase_zpk,...
                          q, delta, N);
% Set Options for Intial Point Optimization
inti_opt_options = optimoptions('ga','UseParallel',true, ...
    'PopulationSize', 5000, ...
    'CrossoverFraction', 0.8, 'MigrationFraction', 0.6,...
    'UseParallel',true, 'InitialPopulationMatrix',init_ystate,...
    'OutputFcn',init_output, 'MaxGenerations', 30,...
    'MaxStallGenerations',10, 'PlotFcn', @gaplotbestf); 
% Run the Genetic Algorithm
disp('---->> Optimization of Initial Values')
if strcmp(option, 'PoleZero')
    init_ystate_opt = ga(@(x)objective(x, 'INIT',...
        q, delta, omega, struture, mag_zpk, phase_zpk,...
        N, Mr, Mc, phi, sd_phi),...
        K+1,[],[],[],[],lb,ub,[],[],inti_opt_options);
else
    init_ystate_opt = ga(@(x)objective(x, 'INIT',...
        q, delta, omega, struture, mag_true, phase_true,...
        N, Mr, Mc, phi, sd_phi),...
        K+1,[],[],[],[],lb,ub,[],[],inti_opt_options);
end
init_kappa = init_ystate_opt(end);       % the optimited kappa
init_y_opt = init_ystate_opt(1:end-1)';  % the optimited values of y
% disp('------>> Optimal solution of Initial Point found by GA solver: ');
% for ii=1:N
%    disp(['-------- R' num2str(ii) '=' num2str(init_y_opt(ii)) ' ohms;', ...
%          'C' num2str(ii) '=' num2str(init_y_opt(N+ii)) ' pF.']);
% end
% disp(['------>> kappa_init=' num2str(init_kappa) '.']);
[~, ~, Xr_zpk, Xc_zpk] = value2code(N, Mr, Mc, init_y_opt, phi, sd_phi);

kappa_zpk = init_kappa;
yr_zpk = init_y_opt(1:N);
yc_zpk = init_y_opt(N+1:end);

if strcmp(option, 'Proposed')
    %% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %%       Optimization with nonlinear contraintment 
    %%             and sparsity regularization
    %% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Create Handle to Custom Output Function
    opt_output = @(a1,a2,a3)optim_plot(a1, a2, a3, omega, struture,...
                              mag_true, phase_true, mag_zpk, phase_zpk,...
                              q, delta, N);
    %  -----------------------------------------------------------------
    %         GA Algorithm --- Barrier function method
    %  ----------------------------------------------------------------- 
    opt_options = optimoptions('ga','UseParallel',true, ...
        'PopulationSize', 2000, 'FunctionTolerance', 1e-6,...
        'UseParallel',true, 'InitialPopulationMatrix',init_ystate_opt,...
        'OutputFcn',opt_output, 'MaxGenerations',30,...
        'MaxStallGenerations',10, 'PlotFcn', @gaplotbestf);
    % Run the Genetic Algorithm
    disp('---->> Sparse Optimization of Component Parameters:')
    ystate_opt = ga(@(x)objective(x, 'ITER', q, delta,...
        omega, struture, mag_true, phase_true, N, Mr, Mc, phi, sd_phi),...
        K+1,[],[],[],[],lb,ub,...
        [],[],opt_options);                     

    % Calculate parameter matrix X 
    kappa = ystate_opt(end); % the optimited kappa
    y_opt = ystate_opt(1:end-1)';  % the optimited values of y

    disp(['------>> Components usage of our method, N=' num2str(N)]);
    disp(['------>> kappa=' num2str(kappa) ';']);
    [yr, yc, Xr, Xc] = value2code(N,Mr,Mc,y_opt,phi,sd_phi);
    for i=1:N
       disp(['-------- R' num2str(i) '=' num2str(yr(i)) ...
             ' ohms,' 'num=' num2str(sum(Xr(i,:))) ';'...
             'C' num2str(i) '=' num2str(yc(i)) ' pF,' ...
             'num=' num2str(sum(Xc(i,:))) '.']);
    end

    % Calculate the maximum discrepancy between \hat{H}(s) and H(s)
    mag_opt = mag_curve(struture, q, delta, omega, ystate_opt, N);
    AvdB = 20*log10(mag_true);
    nAvdB = 20*log10(mag_opt);
    D_sigma = max(abs(nAvdB - AvdB));
    disp(['-------- Maximum discrepancy D_sigma=' num2str(D_sigma) ','...
          'delta=' num2str(delta) '.']);
end
