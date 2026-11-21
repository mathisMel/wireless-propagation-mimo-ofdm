%% Wireless Propagation & MIMO-OFDM — Session 4: Narrowband Local-Area Fading
% Course: ELEC-H422 / H423 (Wireless Communications & Propagation)
%
% Analyzes small-scale spatial fading across a 1 m x 1 m urban-canyon grid:
%   1. Amplitude spatial contours at 3.575 GHz and 3.625 GHz (LOS & NLOS)
%   2. Empirical vs. fitted Rician/Rayleigh PDF & CDF distributions and K-factor
%   3. Spatial autocorrelation and 0.75 coherence distance along X and Y
%   4. Power Angular Spectrum (PAS), RMS angular spread, and spatial array SVDs
% Requires Statistics & Machine Learning Toolbox + Image Processing Toolbox (`mean2`).

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(rootDir, 'functions'));
addpath(fullfile(rootDir, 'src', 'common'));

%% Channel Setup
ch.Ptx = 0.1;
ch.wt = 8;
ch.wb = -8;
ch.rho = 5;
ch.fc = 3.6e9;
ch.B = 100e6;
ch.Q = 1024;
lambda_fc = 3e8 / ch.fc;
dr = lambda_fc / 4;
ch.beta = 2*pi / lambda_fc;
ch.A = ch.Ptx * lambda_fc^2 / (4*pi)^2;
ch.path_params.Npaths = 8;
ch.path_params.max_delay = 500e-9;

%% Amplitude Distribution over Local Area (3.575 vs. 3.625 GHz)
ch.Rx_pos_x = 12:0.1:13;
ch.Rx_pos_y = 6:0.1:7;
carriers = [3.575e9 3.625e9];
[X, Y] = meshgrid(ch.Rx_pos_x, ch.Rx_pos_y);

for i = 1:length(carriers)
    ch.fc = carriers(i);
    lambda_fc = 3e8 / ch.fc;
    ch.beta = 2*pi / lambda_fc;
    ch.A = ch.Ptx * lambda_fc^2 / (4*pi)^2;
    [h_los, h_nlos, ~, ~, ~] = getNarrowBand(ch);

    figure;
    contourf(X, Y, abs(h_los)); colorbar;
    xlabel('x [m]'); ylabel('y [m]');
    title("Amplitude Distribution LOS - " + ch.fc/1e9 + " GHz");

    figure;
    contourf(X, Y, abs(h_nlos)); colorbar;
    xlabel('x [m]'); ylabel('y [m]');
    title("Amplitude Distribution NLOS - " + ch.fc/1e9 + " GHz");
end

%% Fine Grid: K-Factor Extraction + Empirical vs. Fitted PDF / CDF
ch.rho = 1;
ch.fc = 3.6e9;
ch.Rx_pos_x = 12:0.01:13;
ch.Rx_pos_y = 6:0.01:7;
lambda_fc = 3e8 / ch.fc;
ch.beta = 2*pi / lambda_fc;
ch.A = ch.Ptx * lambda_fc^2 / (4*pi)^2;

[h_los, h_nlos, theo_K_los, theo_K_nlos, ~] = getNarrowBand(ch);
P_los = abs(h_los(:));
P_nlos = abs(h_nlos(:));

pdl = fitdist(P_los, 'Rician');
pdnl = fitdist(P_nlos, 'Rician');
pdnl_rayl = fitdist(P_nlos, 'Rayleigh');

num_kl = (pdl.s)^2 / (2*(pdl.sigma^2));
num_knl = (pdnl.s)^2 / (2*(pdnl.sigma^2));

% Report Figs 14-15: Empirical vs. Fitted PDF & CDF
figure;
subplot(1, 2, 1);
histogram(P_los, 40, 'Normalization', 'pdf', 'FaceAlpha', 0.6); hold on;
x_los = linspace(min(P_los), max(P_los), 300);
plot(x_los, pdf(pdl, x_los), 'r-', 'LineWidth', 1.8); hold off; grid on;
xlabel('|h_{LOS}|'); ylabel('PDF'); title('LOS Envelope PDF vs. Rician Fit');
legend('Empirical', sprintf('Rician (K = %.1f dB)', 10*log10(num_kl)), 'Location', 'best');

subplot(1, 2, 2);
histogram(P_nlos, 40, 'Normalization', 'pdf', 'FaceAlpha', 0.6); hold on;
x_nlos = linspace(0, max(P_nlos), 300);
plot(x_nlos, pdf(pdnl, x_nlos), 'r-', 'LineWidth', 1.8);
plot(x_nlos, pdf(pdnl_rayl, x_nlos), 'k--', 'LineWidth', 1.5); hold off; grid on;
xlabel('|h_{NLOS}|'); ylabel('PDF'); title('NLOS Envelope PDF vs. Rician/Rayleigh Fit');
legend('Empirical', sprintf('Rician (K = %.1f dB)', 10*log10(num_knl)), 'Rayleigh', 'Location', 'best');
sgtitle('Small-Scale Fading Amplitude PDFs');

figure;
subplot(1, 2, 1);
[f_los, xi_los] = ecdf(P_los);
plot(xi_los, f_los, 'b-', 'LineWidth', 1.5); hold on;
plot(x_los, cdf(pdl, x_los), 'r--', 'LineWidth', 1.6); hold off; grid on;
xlabel('|h_{LOS}|'); ylabel('CDF'); title('LOS Empirical vs. Rician CDF');
legend('Empirical CDF', 'Rician Fit', 'Location', 'best');

subplot(1, 2, 2);
[f_nlos, xi_nlos] = ecdf(P_nlos);
plot(xi_nlos, f_nlos, 'b-', 'LineWidth', 1.5); hold on;
plot(x_nlos, cdf(pdnl, x_nlos), 'r--', 'LineWidth', 1.6);
plot(x_nlos, cdf(pdnl_rayl, x_nlos), 'k:', 'LineWidth', 1.6); hold off; grid on;
xlabel('|h_{NLOS}|'); ylabel('CDF'); title('NLOS Empirical vs. Rayleigh CDF');
legend('Empirical CDF', 'Rician Fit', 'Rayleigh Fit', 'Location', 'best');
sgtitle('Small-Scale Fading Envelope CDFs');

%% Spatial Autocorrelation & Coherence Distance
nb_real_x = length(ch.Rx_pos_x) - 1;
nb_real_y = length(ch.Rx_pos_y) - 1;
corr_los_x = zeros(1, nb_real_x + 1);
corr_nlos_x = zeros(1, nb_real_x + 1);
corr_los_y = zeros(1, nb_real_y + 1);
corr_nlos_y = zeros(1, nb_real_y + 1);

for i = 0:nb_real_x
    corr_los_x(1, i+1) = mean(h_los(:, i+1:end) .* conj(h_los(:, 1:end-i)), 'all');
    corr_nlos_x(1, i+1) = mean(h_nlos(:, i+1:end) .* conj(h_nlos(:, 1:end-i)), 'all');
end
for i = 0:nb_real_y
    corr_los_y(1, i+1) = mean(h_los(i+1:end, :) .* conj(h_los(1:end-i, :)), 'all');
    corr_nlos_y(1, i+1) = mean(h_nlos(i+1:end, :) .* conj(h_nlos(1:end-i, :)), 'all');
end

corr_los_x = corr_los_x ./ corr_los_x(1, 1);
corr_nlos_x = corr_nlos_x ./ corr_nlos_x(1, 1);
corr_los_y = corr_los_y ./ corr_los_y(1, 1);
corr_nlos_y = corr_nlos_y ./ corr_nlos_y(1, 1);

dx = 0.01;
threshold = 0.75;
idx_los_x = find(abs(corr_los_x) <= threshold, 1);
idx_nlos_x = find(abs(corr_nlos_x) <= threshold, 1);
idx_los_y = find(abs(corr_los_y) <= threshold, 1);
idx_nlos_y = find(abs(corr_nlos_y) <= threshold, 1);

coherence_dist_x_los  = ternary_idx(idx_los_x, nb_real_x, dx);
coherence_dist_x_nlos = ternary_idx(idx_nlos_x, nb_real_x, dx);
coherence_dist_y_los  = ternary_idx(idx_los_y, nb_real_y, dx);
coherence_dist_y_nlos = ternary_idx(idx_nlos_y, nb_real_y, dx);

figure;
plot(0:dx:1, abs(corr_los_x), 'LineWidth', 1.3); hold on;
plot(0:dx:1, abs(corr_nlos_x), 'LineWidth', 1.3); hold off; grid on;
xlabel('\Delta x [m]'); ylabel('|R(\Delta x)|');
title('Spatial Correlation Along X-Axis'); legend('LOS', 'NLOS');

figure;
plot(0:dx:1, abs(corr_los_y), 'LineWidth', 1.3); hold on;
plot(0:dx:1, abs(corr_nlos_y), 'LineWidth', 1.3); hold off; grid on;
xlabel('\Delta y [m]'); ylabel('|R(\Delta y)|');
title('Spatial Correlation Along Y-Axis'); legend('LOS', 'NLOS');

%% Power Angular Spectrum (PAS) & Angular Spread
beta = ch.beta;
theta_list = linspace(0, pi, 180);
z_ind = 1:length(ch.Rx_pos_y);
Z = length(ch.Rx_pos_y);
S_losx = zeros(180, 1); S_nlosx = zeros(180, 1);
S_losy = zeros(180, 1); S_nlosy = zeros(180, 1);
for i = 1:length(theta_list)
    S_losy(i) = corr_los_y * exp(-1j * beta * cos(theta_list(i)) * z_ind ./ Z).';
    S_nlosy(i) = corr_nlos_y * exp(-1j * beta * cos(theta_list(i)) * z_ind ./ Z).';
    S_losx(i) = corr_los_x * exp(-1j * beta * cos(theta_list(i)) * z_ind ./ Z).';
    S_nlosx(i) = corr_nlos_x * exp(-1j * beta * cos(theta_list(i)) * z_ind ./ Z).';
end

figure;
plot(theta_list*180/pi, 10*log10(abs(S_losy))); hold on;
plot(theta_list*180/pi, 10*log10(abs(S_nlosy))); hold off; grid on;
title('Power Angular Spectrum (PAS) — Y Axis');
xlabel('Angle [deg]'); ylabel('Power [dB]'); legend('LOS', 'NLOS');

figure;
plot(theta_list*180/pi, 10*log10(abs(S_losx))); hold on;
plot(theta_list*180/pi, 10*log10(abs(S_nlosx))); hold off; grid on;
title('Power Angular Spectrum (PAS) — X Axis');
xlabel('Angle [deg]'); ylabel('Power [dB]'); legend('LOS', 'NLOS');

fprintf('Session 4 Summary | K_los=%.2f (%.2f dB), K_nlos=%.4f (%.2f dB) | theo_K_los=%.2f dB\n', ...
    num_kl, 10*log10(num_kl), num_knl, 10*log10(num_knl), db(mean2(theo_K_los)));
fprintf('Session 4 Summary | Coherence dist [m] (xLOS, xNLOS, yLOS, yNLOS): %.2f, %.2f, %.2f, %.2f\n', ...
    coherence_dist_x_los, coherence_dist_x_nlos, coherence_dist_y_los, coherence_dist_y_nlos);

function d = ternary_idx(idx, nb_real, dx)
if ~isempty(idx), d = (idx - 1) * dx; else, d = (nb_real - 1) * dx; end
end
