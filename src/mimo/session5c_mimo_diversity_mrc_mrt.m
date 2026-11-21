%% Wireless Propagation & MIMO-OFDM — Session 5c: SIMO MRC & MISO MRT Diversity
% Course: ELEC-H422 / H423 (Wireless Communications & Propagation)
%
% Compares 1xN Receive Maximum Ratio Combining (MRC) and Nx1 Transmit
% Maximum Ratio Transmission (MRT) diversity gains for N = 1..4 antennas over
% Rayleigh fading vs. AWGN QPSK benchmark (Report Fig 21).
%
% Runtime note: Default numSymbols = 2e4, nb_real = 200 runs in ~15-20 s.
% For the full 1000-realization sweep (~3 min), set:
%   numSymbols = 5e4; nb_real = 1000; snr_db = 0:2:40;

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rng(0);
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(rootDir, 'src', 'common'));

%% Signal & Simulation Parameters
numSymbols = 2e4;           % Fast default (report setting: 5e4)
nb_real = 200;              % Fast default (report setting: 1000)
snr_db = 0:2:30;            % Fast default (report setting: 0:2:40)
data_nbps = 2;              % QPSK
min_antenna = 1;
max_antenna = 4;
modulationType = ["qam", "qam"];

BER_MRC  = zeros(max_antenna, length(snr_db), nb_real);
BER_MRT  = zeros(max_antenna, length(snr_db), nb_real);
BER_AWGN = zeros(max_antenna, length(snr_db), nb_real);

h_all = (randn(max_antenna, nb_real) + 1j*randn(max_antenna, nb_real)) / sqrt(2);

%% SIMO MRC & MISO MRT Diversity Sweep
for i = min_antenna:max_antenna
    single_side = 1;
    nb_antenna = i;
    for j = 1:length(snr_db)
        snr = snr_db(j);
        for k = 1:nb_real
            h = h_all(1:nb_antenna, k);
            dataBits = randi([0 1], single_side, data_nbps*numSymbols).';
            x = (mapping_2024(dataBits, data_nbps, modulationType(1), single_side)).';
            power_sig = mean(abs(x).^2, 2);
            power_noise = power_sig / (10^(snr/10));

            n_mrc = sqrt(power_noise/2) .* (randn(nb_antenna, numSymbols) + 1j*randn(nb_antenna, numSymbols));
            n_mrt = sqrt(power_noise/2) .* (randn(single_side, numSymbols) + 1j*randn(single_side, numSymbols));

            y_awgn = x + n_mrc;

            % 1. Receive Diversity: Maximum Ratio Combining (MRC)
            y_mrc = h * x + n_mrc;
            x_MRC = h' * y_mrc;

            % 2. Transmit Diversity: Maximum Ratio Transmission (MRT)
            h_los_h = h.' * conj(h);
            y_MRT = (h_los_h / sqrt(h_los_h)) * x + n_mrt;

            dataBits_MRC  = demapping(x_MRC.', data_nbps, modulationType(1), single_side);
            dataBits_MRT  = demapping(y_MRT.', data_nbps, modulationType(1), single_side);
            dataBits_AWGN = demapping(y_awgn.', data_nbps, modulationType(1), nb_antenna);
            dataBits_compare = repmat(dataBits(:), 1, nb_antenna);

            BER_MRC(i, j, k)  = norm(dataBits(:) - dataBits_MRC(:), 1) / numel(dataBits_MRC);
            BER_MRT(i, j, k)  = norm(dataBits(:) - dataBits_MRT(:), 1) / numel(dataBits_MRT);
            BER_AWGN(i, j, k) = norm(dataBits_compare(:) - dataBits_AWGN(:), 1) / numel(dataBits_AWGN);
        end
    end
end

BER_MRC  = mean(BER_MRC, 3);
BER_MRT  = mean(BER_MRT, 3);
BER_AWGN = mean(mean(BER_AWGN, 3), 1);

%% Plot Diversity Gains (Report Fig 21)
figure;
semilogy(snr_db, BER_MRC, 'LineWidth', 1.2); hold on;
semilogy(snr_db, BER_AWGN, 'k-o', 'LineWidth', 1.2);
semilogy(snr_db, BER_MRT, '--', 'LineWidth', 1.2); hold off;
xlabel('SNR [dB]'); ylabel('BER');
title('BER vs. SNR for SIMO MRC and MISO MRT Diversity');
legend('MRC Rx=1', 'MRC Rx=2', 'MRC Rx=3', 'MRC Rx=4', 'AWGN', ...
       'MRT Tx=1', 'MRT Tx=2', 'MRT Tx=3', 'MRT Tx=4', 'Location', 'southwest');
grid on;
set(gca, 'YScale', 'log');
ylim([1e-5 1]);

fprintf('Session 5c @ 14 dB | MRC(1..4): %s | AWGN: %.3e\n', ...
    sprintf('%.3e ', BER_MRC(:, snr_db==14)), BER_AWGN(snr_db==14));
