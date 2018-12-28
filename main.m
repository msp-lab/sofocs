%% Find the Optimal Component Values for a Chaos Circuit

%  Copyright (c) 2018, Kunpeng Wang.
%  Email: wkphnzk@163.com

%% WARNING!
% 1. DO NOT DELETE FILE: stdrc.mat
% 2. IF YOU WANT TO RECALCULATE CIRCUIT PARAMETERS, DELETE THE
%    CORRESPONDING FILES IN DATA DIRECTORY
% 3. THEN, CHANGE THE PARAMETERS IN "(2) Set parameters"

%% Step 1: Load Parameters
close all; clear;
data_dir = 'data/';
if ~exist(data_dir,'dir'), mkdir(data_dir); end

disp('Step 1: Initial Experiment Parameters')
% (1) Load in standard component values
% Resistor  always choose E96 1% tolerance
% Capicitor always choose E12 10% tolerance
disp('++ Load in standard component values')
load('stdrc.mat', 'rstd', 'cstd', 'sd_rstd', 'Mr', 'Mc',...
                  'sd_cstd', 'omega', 'phi', 'sd_phi');
         
%% (2) Set parameters
disp('++ Set simulation parameters')
q = 0.75;    % fractional order
delta = 1;  % dB, The maximum discrepancy
structure = 'Chain'; % Fractance structure {'Chain', 'Tree',  'Ladder'}
N = 4; % system order of franctance

MaxStall = 10;
MaxUC = 0.2; 
UCItNum = 2000;
result_fn = [data_dir structure '_q' num2str(q) '_N' num2str(N) '_' ...
             num2str(delta) 'dB'];

%% (3) Start processing
disp('Start processing')
% get the system order of transfer function
[mag_true, phase_true] = mag_curve('Desired', q, delta, omega);
[mag_zpk, phase_zpk] = mag_curve('ZPK', q, delta, omega);

                                   
%% loop no more than it_max times
if ~exist([result_fn '.mat'], 'file')
    iter = 1;
    state_calc = 0;
    UC = 1.;
    [~, ~, N_zpk] = mag_curve('ZPK', q, delta, omega);
    [Xr_zpk, Xc_zpk, yr_zpk, yc_zpk, kappa_zpk] = ...
        sofcs(N_zpk, omega, structure, mag_true, phase_true, ...
              mag_zpk, phase_zpk, Mr, Mc, q, delta, phi, sd_phi, ...
              'Zigzag');
    [UC_zpk, ~] = uncertainty_eval(N_zpk, structure, kappa_zpk, q, delta,...
                                   Xr_zpk, Xc_zpk, UCItNum, data_dir);
    while iter <= MaxStall
        %% Start sparse optimization
        disp(['++ The maximum discrepancy delta=' num2str(delta) 'dB'])
        disp(['++-- Fractance strcture: ' structure])
        disp(['++-- Franctance order N=' num2str(N) ',q=' num2str(q)])
        [~, ~, ~, ~, ~, kappa, Xr, Xc, yr, yc, D_sigma] = ...
            sofcs(N, omega, structure, mag_true, phase_true, ...
                  mag_zpk, phase_zpk, Mr, Mc, q, delta, ...
                  phi, sd_phi, 'Proposed');
              
        disp(['---- @@ iter=' num2str(iter)])
        iter = iter + 1;        
        %% Save the minimum requirement system order N_rc
        if D_sigma <= delta                                      
            [UC, ~] = uncertainty_eval(N, structure, kappa, q, delta,...
                                       Xr, Xc, UCItNum, data_dir);
            if UC < 0.2    
                state_calc = 1;
                disp(['Success! iter=' num2str(iter)]);
                % Save results     
                save([result_fn '.mat'], 'q', 'N', 'structure', 'delta',...
                    'kappa', 'Xr', 'Xc', 'yr', 'yc', 'D_sigma', ...
                    'Xr_zpk', 'Xc_zpk', 'yr_zpk', 'yc_zpk', 'kappa_zpk');
                break;      
            end
        else % D_sigma > delta
            disp(['++-- D_sigma > delta, termination!' ...
                  num2str(D_sigma) '>' num2str(delta)])
        end %-- if D_sigma
    end
else
    load([result_fn '.mat'])
    state_calc = 1;
end

disp(['++-- @@ ', structure, '_', num2str(delta), 'dB finished.'])
if ~state_calc
    disp(['Calculation Failure! UC=' num2str(UC) ...
          'You can increase the system order of fractance (N), ' ...
          'and then try again.']);
else
    if exist([result_fn '.txt'], 'file')
        delete([result_fn '.txt'])
    end
    fwhd = fopen([result_fn '.txt'], 'w+');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Zigzag method
    info = ['-->> Components usage of Zigzag method, N=' num2str(N_zpk) '\n'...
            '++-- The maximum discrepancy delta=' num2str(delta) 'dB\n' ...
            '++-- Fractance strcture: ' structure '\n'...
            '++-- Franctance order N=' num2str(N_zpk) ',q=' num2str(q) '\n'];
    fprintf(info);
    fprintf(fwhd, info);
    for i=1:N_zpk
        result = ['-------- R' num2str(i) '=' sprintf('%9d', yr_zpk(i)') ...
            '\tSeries   : ' ...
             sprintf('%9d', rstd(Xr_zpk(i,:)>0)')  ' ohms, ' ...
             'num=' num2str(sum(Xr_zpk(i,:))) ';\n'];
       fprintf(result);
       fprintf(fwhd, [result '\n']);
       result = ['-------- C' num2str(i) '=' sprintf('%9d', yc_zpk(i)) ...
           '\tParallel : ' ...
             sprintf('%9d', cstd(Xc_zpk(i,:)>0)')  ' pF,   ' ...
             'num=' num2str(sum(Xc_zpk(i,:))) '.\n'];
       fprintf(result);
       fprintf(fwhd, [result '\n']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % our method
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
    info = ['-->> Components usage of our method, N=' num2str(N) '\n'...
            '++-- The maximum discrepancy delta=' num2str(delta) 'dB\n' ...
            '++-- Fractance strcture: ' structure '\n'...
            '++-- Franctance order N=' num2str(N) ',q=' num2str(q) '\n'];
    fprintf(info);
    fprintf(fwhd, '+++++++++++++++++++++++++++++++++++++++++++++++++++\n');
    fprintf(fwhd, info);
    for i=1:N
        result = ['-------- R' num2str(i) '=' sprintf('%9d', yr(i)') ...
            '\tSeries   : ' ...
             sprintf('%9d', rstd(Xr(i,:)>0)')  ' ohms, ' ...
             'num=' num2str(sum(Xr(i,:))) ';\n'];
       fprintf(result);
       fprintf(fwhd, [result '\n']);
       result = ['-------- C' num2str(i) '=' sprintf('%9d', yc(i)) ...
           '\tParallel : ' ...
             sprintf('%9d', cstd(Xc(i,:)>0)')  ' pF,   ' ...
             'num=' num2str(sum(Xc(i,:))) '.\n'];
       fprintf(result);
       fprintf(fwhd, [result '\n']);
    end
    fclose(fwhd);
end






