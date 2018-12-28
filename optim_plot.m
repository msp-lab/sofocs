function [state,options,optchanged] = optim_plot(options, state, flag,...
    omega, struture, mag_true, phase_true, mag_cmp, phase_cmp,...
    q, delta, N)
%% Plot the Result of Optimal Values
%% Input:
%  omega: vector of frequency sample points
%  struture : The structure of fractance, Desired/Chain/Tree/Ladder/ZPK
%  mag_true : vector of desired Magnitude Response Curve
%  phase_true: vector of desired phase of Magnitude Response Curve
%  mag_cmp :  Magnitude Response Curve of zpk method
%  phase_cmp: phase of Magnitude Response Curve of zpk method
%  q : the order of fractional-order chaos system
%  delta : dB, the maximum discrepancy between desired H(s) and \hat{H}(s)
%  N: system order of franctance 
%
% Copyright (c) 2018, Kunpeng Wang.
% Email: wkphnzk@163.com
%

optchanged = false;

fig_name = 'MagResCurve'; 
plot_name = 'SOFCS ';

switch flag
    case 'init'       
        fig = findobj(0,'Tag',fig_name);
        if isempty(fig)
            figure('Position',[400 100 800 600], 'Name', plot_name,...
                   'NumberTitle', 'off', 'Tag', fig_name);
        else
            figure(fig(end));
        end
        
    case {'iter', 'done'}
        fig = findobj(0,'Tag',fig_name);
        if isempty(fig)
            figure('Position',[400 100 800 600], 'Name', plot_name,...
                   'NumberTitle', 'off', 'Tag', fig_name);
        else
            figure(fig(end));
        end

        [~,loc] = min(state.Score); % Find location of best
        [best_mag, best_phase] = mag_curve(struture,q,delta,omega,...
                     state.Population(loc,:),N);
                        
        subplot(2,1,1); 
        semilogx(omega,20*log10(mag_true),'-r'); hold on; grid on;
        semilogx(omega,20*log10(mag_cmp),'-.m');  % comparision curve
        semilogx(omega,20*log10(best_mag),'--b'); % plot the GA curve
        legend('Ideal','Comparison',[plot_name 'Solution']);
        title([plot_name 'Optimization: Magnitude Response Curve']); 
        ylabel('|H(j\omega)|'); 
        hold off;
        subplot(2,1,2); 
        semilogx(omega,phase_true,'-r'); hold on; grid on;
        semilogx(omega,phase_cmp,'-.m');
        semilogx(omega,best_phase,'--b'); % plot the GA curve
        xlabel('\omega (rad/sec)'); ylabel('\angleH(j\omega) (\circ)'); 
        legend('Ideal','Comparison',[plot_name 'Solution']);
        title([plot_name 'Optimization: Phase Response Curve']);
        drawnow;
        hold off;
end
