function [h] = plotEMCDfitCum(fstruct, varargin)
%PLOTEMCDFIT Uses the output structure of getMLMSfromModel.m to generate a
%plot.  Returns a handle

%% Set default values
xlim = [660 860];

%% Get input pair values
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for idx = 1:2:(length(varargin)-1)
        switch lower(varargin{idx})
            case 'xlim'
                xlim = varargin{idx + 1};
        end
    end
end

%% Begin

%   Get figure size
monitors = get(0, 'MonitorPositions');
width = monitors(1, 3) .* .95;
ratio = 0.1;
% pos = [0 0 width width/ratio];
pos = [0 0 1200 800];


% Center figure
h.fig = figure('color', 'w', 'position', pos);
movegui(h.fig,'onscreen')

% f = fstruct.f;
e = fstruct.e;
% plus = fstruct.preplus + f(1) .* e .^ f(2);
% minus = fstruct.preminus + f(3) .* e .^ f(4);

%% Raw data plot
h.ax(1) = subplot(2, 2, 1);
h.pl(1) = plot(e, fstruct.plus, 'DisplayName', 'Plus', 'LineWidth', 2);
hold on;
h.pl(2) = plot(e, fstruct.minus, 'DisplayName', 'Minus', 'LineWidth', 2);
h.pl(3) = plot(e, fstruct.BGplus, 'r--', 'DisplayName', 'BG Plus');
h.pl(4) = plot(e, fstruct.BGminus, 'r--', 'DisplayName', 'BG Minus');
legend;
set(gca, 'xlim', xlim, 'FontSize', 14);
h.ax(1).Tag = 'Raw Plot';
xlabel('Energy [eV]');
ylabel('Intensity [counts]');
grid on;
title('Raw spectra');
% linkdata on;

%% EMCD signal plot
h.ax(2) = subplot(2, 2, 2);
h.pl(5) = plot(e, fstruct.preplus, 'DisplayName', 'Plus BGsub', 'LineWidth', 1.5);
hold on;
h.pl(6) = plot(e, fstruct.preminus, 'DisplayName', 'Minus BGsub', 'LineWidth', 1.5);

h.pl(7) = plot(e, fstruct.EMCD, 'k', 'DisplayName', 'EMCD signal', 'LineWidth', 2);
h.pl(8) = plot(e, fstruct.EMCDsig, 'r', 'DisplayName', 'EMCD fit');
legend;
set(gca, 'xlim', xlim, 'FontSize', 14);
h.ax(2).Tag = 'Norm Plot';
xlabel('Energy [eV]');
ylabel('Intensity [normalized]');
grid on;
% title('Normalized and BG subtracted');
title(sprintf('L_3 = %2.2f %%, L_2 = %2.2f %%', abs(fstruct.f(13)*100), abs(fstruct.f(19)*100)));

%% Plot residuals
h.ax(3) = subplot(2, 2, 3);
% set(h.ax(3), 'Tag', 'EMCD Plot');
h.pl(9) = plot(e, fstruct.EMCD, 'k', 'DisplayName', 'EMCD signal', 'LineWidth', 2);
hold on;
h.pl(10) = plot(e, fstruct.EMCDsig, 'r', 'DisplayName', 'EMCD fit');
legend;
set(gca, 'xlim', xlim, 'FontSize', 14);
h.ax(3).Tag = 'EMCD Plot';
xlabel('Energy [eV]');
ylabel('Intensity [percent of L_3]');
grid on;
% title(sprintf('Residuals, m_L/m_S = %0.4f', fstruct.mLmS));
title(sprintf('SNR L_3 = %2.2f dB, L_2 = %2.2f dB', fstruct.SNRL3, fstruct.SNRL2));

%% Plot cumulative sum
h.ax(4) = subplot(2, 2, 4);
% set(h.ax(3), 'Tag', 'EMCD Plot');
h.pl(11) = plot(e(makeWin(e, xlim)), ...
    cumsum(fstruct.EMCD(makeWin(e, xlim))), 'k', ...
    'DisplayName', 'EMCD cumsum', 'LineWidth', 2);
hold on;
h.pl(12) = plot(e(makeWin(e, xlim)), ...
    cumsum(fstruct.EMCDsig(makeWin(e, xlim))), 'r', ...
    'DisplayName', 'EMCD fit cumsum');
% legend;
set(gca, 'xlim', xlim, 'FontSize', 14);
h.ax(4).Tag = 'Cumsum Plot';
xlabel('Energy [eV]');
ylabel('Intensity');
grid on;
title(sprintf('Cumulative Sum, m_L/m_S = %0.4f', fstruct.mLmS));

end

