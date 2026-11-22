%% Wireless Propagation & MIMO-OFDM — Session 5a: Narrowband MIMO Singular Values
% Course: ELEC-H422 / H423 (Wireless Communications & Propagation)
%
% Inspects 1x4 SIMO / MIMO channel correlation matrix singular values in the
% urban street-canyon model along both Y-oriented and X-oriented arrays (LOS/NLOS).

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(rootDir, 'functions'));
addpath(fullfile(rootDir, 'src', 'common'));

%% Common Channel Parameters
ch.Ptx = 0.1;
ch.wt = 8;
ch.wb = -8;
ch.fc = 3.6e9;
ch.B = 100e6;
ch.Q = 1024;
lambda_fc = 3e8 / ch.fc;
dr = lambda_fc / 4;
ch.beta = 2*pi / lambda_fc;
ch.A = ch.Ptx * lambda_fc^2 / (4*pi)^2;
ch.Rx_pos_x = 12;
ch.Rx_pos_y = 6;
ch.Ntx = 1;
ch.Nrx = 4;
ch.AntSpacing = dr;

%% 1. Y-Oriented Receive Array (rho = 10, Npaths = 8)
ch.rho = 10;
ch.path_params.Npaths = 8;
ch.path_params.max_delay = 500e-9;
ch.RxArrayOrient = 'y';
[h_los_y, h_nlos_y, ~] = getNarrowBandMIMO(ch);
h_los_y = h_los_y ./ max(abs(h_los_y));
h_nlos_y = h_nlos_y ./ max(abs(h_nlos_y));
[~, S1, ~] = svd(h_los_y * h_los_y');
[~, S2, ~] = svd(h_nlos_y * h_nlos_y');

%% 2. X-Oriented Receive Array (rho = 3, Npaths = 5)
ch.rho = 3;
ch.path_params.Npaths = 5;
ch.RxArrayOrient = 'x';
[h_los_x, h_nlos_x, ~] = getNarrowBandMIMO(ch);
h_los_x = h_los_x ./ max(abs(h_los_x));
h_nlos_x = h_nlos_x ./ max(abs(h_nlos_x));
[~, S_los, ~] = svd(h_los_x * h_los_x');
[~, S_nlos, ~] = svd(h_nlos_x * h_nlos_x');

sv_mat = [diag(S1), diag(S2), diag(S_los), diag(S_nlos)];
figure;
bar(sv_mat);
set(gca, 'YScale', 'log');
grid on;
xlabel('Singular Value Index'); ylabel('Singular Value Magnitude (log scale)');
title('Narrowband 1\times4 Channel Correlation Matrix Singular Values');
legend('Y-array LOS', 'Y-array NLOS', 'X-array LOS', 'X-array NLOS', 'Location', 'northeast');

fprintf('Session 5a Singular Values:\n');
fprintf('  Y-array LOS  : %s\n', sprintf('%.3e ', diag(S1)));
fprintf('  Y-array NLOS : %s\n', sprintf('%.3e ', diag(S2)));
fprintf('  X-array LOS  : %s\n', sprintf('%.3e ', diag(S_los)));
fprintf('  X-array NLOS : %s\n', sprintf('%.3e ', diag(S_nlos)));

%% Export Figures (450 DPI PNG + Vector PDF)
figRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'docs', 'figures');
pngDir = fullfile(figRoot, 'png');
pdfDir = fullfile(figRoot, 'pdf');
if ~exist(pngDir, 'dir'), mkdir(pngDir); end
if ~exist(pdfDir, 'dir'), mkdir(pdfDir); end
figs = findall(0, 'Type', 'figure');
if ~isempty(figs)
    [~, ord] = sort([figs.Number]);
    figs = figs(ord);
    for fi = 1:numel(figs)
        exportgraphics(figs(fi), fullfile(pngDir, sprintf('session5a_svd_fig%02d.png', fi)), 'Resolution', 450, 'BackgroundColor', 'white');
        exportgraphics(figs(fi), fullfile(pdfDir, sprintf('session5a_svd_fig%02d.pdf', fi)), 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
