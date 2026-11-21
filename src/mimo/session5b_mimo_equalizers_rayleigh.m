%% Wireless Propagation & MIMO-OFDM — Session 5b: Rayleigh MIMO Equalizers (MF, ZF, MMSE, ML)
% Course: ELEC-H422 / H423 (Wireless Communications & Propagation)
%
% Compares Matched Filter (MF), Zero-Forcing (ZF), Minimum Mean Square Error (MMSE),
% and Maximum Likelihood (ML) receivers for i x i Rayleigh MIMO (i = 1..4, QPSK).
%
% Runtime note: Default numSymbols = 2e3, nb_real = 10 runs in ~5-8 s.
% For the full high-precision report sweep (slow, ~10 min), set:
%   numSymbols = 5e4; snr_db = 0:2:40;

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rng(0);
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(rootDir, 'src', 'common'));

%% Signal & Simulation Parameters
numSymbols = 2e3;           % Fast default (report setting: 5e4)
nb_real = 10;
snr_db = 0:5:35;            % Fast default (report setting: 0:2:40)
data_nbps = 2;              % QPSK
qpsk_symbols = [1+1j, 1-1j, -1+1j, -1-1j] / sqrt(2);
MIMO_min = 1;
MIMO_max = 4;
modulationType = ["qam", "qam"];

BER_MF   = zeros(MIMO_max, length(snr_db), nb_real);
BER_MMSE = zeros(MIMO_max, length(snr_db), nb_real);
BER_ZF   = zeros(MIMO_max, length(snr_db), nb_real);
BER_ML   = zeros(MIMO_max, length(snr_db), nb_real);

%% MIMO Equalization Loop Across Array Sizes 1x1..4x4
for i = MIMO_min:MIMO_max
    cellArray = repmat({qpsk_symbols}, 1, i);
    gridArray = cell(1, i);
    [gridArray{1:i}] = ndgrid(cellArray{:});
    combos = zeros(i, numel(gridArray{1}));
    for b = 1:i
        combos(b, :) = gridArray{b}(:).';
    end

    H_many = (randn(i, i, nb_real) + 1j*randn(i, i, nb_real)) / sqrt(2);
    avg_energy = mean(abs(H_many), 'all');
    H_many = H_many / abs(avg_energy);

    for j = 1:length(snr_db)
        snr = snr_db(j);
        for k = 1:nb_real
            dataBits = randi([0 1], i, data_nbps*numSymbols).';
            x = (mapping_2024(dataBits, data_nbps, modulationType(1), i)).';
            power_sig = mean(abs(x).^2, 2);
            power_noise = power_sig / (10^(snr/10));
            n = sqrt(power_noise/2) .* (randn(i, numSymbols) + 1j*randn(i, numSymbols));

            y = H_many(:, :, k) * x + n;

            % 1. Matched Filter (MF)
            F_mf = diag(1 ./ sum(abs(H_many(:, :, k)).^2, 1)) * H_many(:, :, k)';
            x_mf = F_mf * y;

            % 2. Zero-Forcing (ZF)
            F_zf = pinv(H_many(:, :, k));
            x_zf = F_zf * y;

            % 3. Maximum Likelihood (ML)
            H_all = H_many(:, :, k) * combos;
            x_ml_hat = zeros(i, numSymbols);
            for m = 1:numSymbols
                diff = H_all - y(:, m);
                dists = sum(abs(diff).^2, 1);
                [~, min_idx] = min(dists);
                x_ml_hat(:, m) = combos(:, min_idx);
            end

            % 4. Minimum Mean Square Error (MMSE)
            I_tx = eye(i);
            F_mmse = (H_many(:, :, k)' * H_many(:, :, k) + power_noise .* I_tx) \ H_many(:, :, k)';
            x_mmse = F_mmse * y;

            dataBits_MF   = demapping(x_mf.', data_nbps, modulationType(1), i);
            dataBits_MMSE = demapping(x_mmse.', data_nbps, modulationType(1), i);
            dataBits_ZF   = demapping(x_zf.', data_nbps, modulationType(1), i);
            dataBits_ML   = demapping(x_ml_hat.', data_nbps, modulationType(1), i);

            BER_MF(i, j, k)   = norm(dataBits(:) - dataBits_MF(:), 1) / numel(dataBits_MF);
            BER_MMSE(i, j, k) = norm(dataBits(:) - dataBits_MMSE(:), 1) / numel(dataBits_MMSE);
            BER_ZF(i, j, k)   = norm(dataBits(:) - dataBits_ZF(:), 1) / numel(dataBits_ZF);
            BER_ML(i, j, k)   = norm(dataBits(:) - dataBits_ML(:), 1) / numel(dataBits_ML);
        end
    end
end

BER_MF   = mean(BER_MF, 3);
BER_MMSE = mean(BER_MMSE, 3);
BER_ZF   = mean(BER_ZF, 3);
BER_ML   = mean(BER_ML, 3);

%% Plots (Report Fig 20)
figure;
for i = MIMO_min:MIMO_max
    subplot(2, 2, i);
    semilogy(snr_db, BER_MF(i, :), '-^'); hold on;
    semilogy(snr_db, BER_MMSE(i, :), '-s');
    semilogy(snr_db, BER_ZF(i, :), '-o');
    semilogy(snr_db, BER_ML(i, :), '-*'); hold off;
    xlabel('SNR [dB]'); ylabel('BER');
    title(sprintf('MIMO N_{tx} = %d, N_{rx} = %d', i, i));
    legend('MF', 'MMSE', 'ZF', 'ML', 'Location', 'southwest');
    grid on;
end
sgtitle('Rayleigh MIMO Equalizer Comparison (QPSK)');

fprintf('Session 5b 4x4 MIMO BER @ 20 dB (MF, ZF, MMSE, ML): %.3e, %.3e, %.3e, %.3e\n', ...
    BER_MF(4, snr_db==20), BER_ZF(4, snr_db==20), BER_MMSE(4, snr_db==20), BER_ML(4, snr_db==20));
