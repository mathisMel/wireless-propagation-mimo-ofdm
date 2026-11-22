%% Wireless Propagation & MIMO-OFDM — Sessions 2 & 3: Complete OFDM Transceiver
% Course: ELEC-H422 / H423 (Wireless Communications & Propagation)
%
% Simulates an end-to-end 1024-subcarrier 16-QAM OFDM link across an urban
% street-canyon wideband channel (LOS & NLOS, CP = 64 samples):
%   - Schmidl-Cox repetitive training preamble with random frame arrival offset
%   - Preamble channel estimation and CTF MSE vs. SNR
%   - Zero-forcing frequency-domain equalization vs. unequalized BER
% Derived from the verified March 2025 project archive.

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rng(0);
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(rootDir, 'functions'));
addpath(fullfile(rootDir, 'src', 'common'));

%% Channel Parameters
ch.Ptx = 0.1;
ch.wt = 8;
ch.wb = -8;
ch.rho = 10;
ch.fc = 3.6e9;
ch.B = 100e6;
ch.Q = 1024;
lambda_fc = 3e8 / ch.fc;
ch.beta = 2*pi / lambda_fc;
ch.A = ch.Ptx * lambda_fc^2 / (4*pi)^2;
ch.path_params.Npaths = 5;
ch.Rx_pos_x = 12;
ch.Rx_pos_y = 6;
Ts = 1 / ch.B;

%% OFDM Transceiver Parameters
num_0_sc = 20;
pilot_percentage = 0.4;
data_nbps = 4;              % 16-QAM data
pilot_nbps = 1;             % BPSK pilots
num_ofdm_blocks = 100;
size_cp = 64;
nb_source = 1;
numNon0SC = 1/2 * (ch.Q - 2*num_0_sc - 1 - 3);
modulationType = ["qam", "bpsk"];

%% Wideband CIR & CTF Extraction
[H_los, H_nlos, ~] = getWideBand_open(ch);
h_los = ifft(H_los, [], 2);
h_nlos = ifft(H_nlos, [], 2);

crop_factor = 300;
tau = (0:crop_factor-1) * Ts;
PDP_los = abs(h_los(1, 1:crop_factor)).^2;
P_los = sum(PDP_los);
tau_mean_los = (1 / P_los) * sum(PDP_los .* tau);
sigma_tau_los = sqrt((1 / P_los) * sum(PDP_los .* tau.^2) - tau_mean_los^2);

size_channel = ceil(sigma_tau_los / Ts);
toa_error_estimate = size_cp - size_channel;

h_los = h_los(1, 1:size_channel) ./ max(abs(h_los(1, 1:size_channel)));
h_nlos = h_nlos(1, 1:size_channel) ./ max(abs(h_nlos(1, 1:size_channel)));
H_comp_los = fft(h_los, ch.Q);
H_comp_nlos = fft(h_nlos, ch.Q);

figure;
subplot(2, 2, 1); stem(abs(h_los)); title('Normalized CIR — LOS'); xlabel('Tap'); ylabel('|h|'); grid on;
subplot(2, 2, 2); stem(abs(h_nlos)); title('Normalized CIR — NLOS'); xlabel('Tap'); ylabel('|h|'); grid on;
subplot(2, 2, 3); plot(abs(H_comp_los)); title('CTF Magnitude — LOS'); xlabel('Subcarrier'); ylabel('|H|'); grid on;
subplot(2, 2, 4); plot(abs(H_comp_nlos)); title('CTF Magnitude — NLOS'); xlabel('Subcarrier'); ylabel('|H|'); grid on;
sgtitle('Urban Street-Canyon Wideband Channel Used for OFDM Simulation');

%% SNR Sweep
% Note: Increase nb_real (e.g. 50-100) for smoother low-SNR curves.
nb_real = 3;
snr_db = -10:5:30;

BER_notEqualizedLos = zeros(length(snr_db), nb_real);
BER_EqualizedLos    = zeros(length(snr_db), nb_real);
BER_notEqualizedNlos = zeros(length(snr_db), nb_real);
BER_EqualizedNlos   = zeros(length(snr_db), nb_real);
toa_error_los       = zeros(length(snr_db), nb_real);
toa_error_nlos      = zeros(length(snr_db), nb_real);
mse_Los             = zeros(length(snr_db), nb_real);
mse_Nlos            = zeros(length(snr_db), nb_real);

window_size = 2 * ch.Q;
preamble_cp_size = 2 * size_cp;

for i = 1:length(snr_db)
    for j = 1:nb_real
        dataBits = randi([0 1], 1, data_nbps*num_ofdm_blocks*ch.Q).';
        pilotBits = randi([0 1], 1, num_ofdm_blocks*pilot_nbps).';
        mappedPilots = mapping_2024(pilotBits, pilot_nbps, modulationType(2), nb_source);
        mappedData = mapping_2024(dataBits, data_nbps, modulationType(1), nb_source);

        [OFDMBlocksTx, dataIndices, pilotIndices, ~] = CreateOFDMBlocks_ELEC422_423( ...
            mappedData, mappedPilots, pilot_percentage, ch.Q, num_0_sc, num_ofdm_blocks, nb_source);

        OFDMBlocksTx = [zeros(num_ofdm_blocks, num_0_sc, nb_source) OFDMBlocksTx(:, 1:end/2, :) ...
                        zeros(num_ofdm_blocks, 3, nb_source) OFDMBlocksTx(:, end/2+1:end, :) ...
                        zeros(num_ofdm_blocks, num_0_sc+1, nb_source)];

        preamble = [zeros(1, num_0_sc, nb_source) sign(randn(1, numNon0SC)) ...
                    zeros(1, 3, nb_source) sign(randn(1, numNon0SC)) zeros(1, num_0_sc+1, nb_source)];
        preamble_time = ifft(preamble, [], 2);
        preamble_time_cp = [preamble_time(:, end-2*size_cp+1:end) preamble_time preamble_time];

        timeBlocksTx = ifft(OFDMBlocksTx, [], 2);
        timeBlocksTx_cp = [timeBlocksTx(:, end-size_cp+1:end, :) timeBlocksTx];
        timeBlocksTx_cp = reshape(permute(timeBlocksTx_cp, [2 1 3]), 1, [], nb_source);
        timeSignalTx = permute(timeBlocksTx_cp, [3 2 1]);
        timeSignalTx_withPreamble = [preamble_time_cp timeSignalTx];
        size_Signaltx = length(timeSignalTx_withPreamble);

        toa_uncertainty = randi(10000);
        timeSignalTx_Final = [zeros(1, toa_uncertainty) timeSignalTx_withPreamble];

        showPlateau = (i == length(snr_db) && j == 1);
        [BER_EqualizedLos(i, j), BER_notEqualizedLos(i, j), mse_Los(i, j), toa_error_los(i, j)] = ...
            BerSimulation_autoCorr(toa_uncertainty, size_Signaltx, timeSignalTx_Final, h_los, snr_db(i), ...
            size_cp, H_comp_los, dataIndices, num_0_sc, data_nbps, modulationType, dataBits, mappedData, ...
            mappedPilots, ch.Q, preamble, pilotIndices, window_size, preamble_cp_size, toa_error_estimate, showPlateau);

        [BER_EqualizedNlos(i, j), BER_notEqualizedNlos(i, j), mse_Nlos(i, j), toa_error_nlos(i, j)] = ...
            BerSimulation_autoCorr(toa_uncertainty, size_Signaltx, timeSignalTx_Final, h_nlos, snr_db(i), ...
            size_cp, H_comp_nlos, dataIndices, num_0_sc, data_nbps, modulationType, dataBits, mappedData, ...
            mappedPilots, ch.Q, preamble, pilotIndices, window_size, preamble_cp_size, toa_error_estimate, false);
    end
end

BER_avgEqualized_los     = mean(BER_EqualizedLos, 2);
BER_avgNotEqualized_los  = mean(BER_notEqualizedLos, 2);
mse_avg_los              = mean(mse_Los, 2);
BER_avgEqualized_nlos    = mean(BER_EqualizedNlos, 2);
BER_avgNotEqualized_nlos = mean(BER_notEqualizedNlos, 2);
mse_avg_nlos             = mean(mse_Nlos, 2);

%% Plots (Report Fig 11)
figure;
subplot(2, 3, 1); semilogy(snr_db, mse_avg_los, '-o'); grid on;
title('Channel Estimation MSE — LOS'); xlabel('SNR [dB]'); ylabel('MSE');

subplot(2, 3, 2); semilogy(snr_db, BER_avgNotEqualized_los, '-o'); grid on;
title('BER Without Equalization — LOS'); xlabel('SNR [dB]'); ylabel('BER');

subplot(2, 3, 3); semilogy(snr_db, BER_avgEqualized_los, '-o'); grid on;
title('BER With ZF Equalization — LOS'); xlabel('SNR [dB]'); ylabel('BER');

subplot(2, 3, 4); semilogy(snr_db, mse_avg_nlos, '-o'); grid on;
title('Channel Estimation MSE — NLOS'); xlabel('SNR [dB]'); ylabel('MSE');

subplot(2, 3, 5); semilogy(snr_db, BER_avgNotEqualized_nlos, '-o'); grid on;
title('BER Without Equalization — NLOS'); xlabel('SNR [dB]'); ylabel('BER');

subplot(2, 3, 6); semilogy(snr_db, BER_avgEqualized_nlos, '-o'); grid on;
title('BER With ZF Equalization — NLOS'); xlabel('SNR [dB]'); ylabel('BER');
sgtitle('16-QAM OFDM Performance Over Urban-Canyon Channel (CP = 64)');

fprintf('Session 2/3 Summary | LOS Eq BER (-10..30 dB): %s\n', sprintf('%.3e ', BER_avgEqualized_los));
fprintf('Session 2/3 Summary | LOS CTF MSE (-10..30 dB): %s\n', sprintf('%.3e ', mse_avg_los));

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
        exportgraphics(figs(fi), fullfile(pngDir, sprintf('session23_ofdm_fig%02d.png', fi)), 'Resolution', 450, 'BackgroundColor', 'white');
        exportgraphics(figs(fi), fullfile(pdfDir, sprintf('session23_ofdm_fig%02d.pdf', fi)), 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
