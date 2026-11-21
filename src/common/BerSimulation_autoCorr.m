function [BER_Equalized, BER_notEqualized, mse, ToA_error] = BerSimulation_autoCorr(true_ToA, size_Signaltx, timeSignalTx, h, snr, size_cp, H, dataIndices, num_0_sc, data_nbps, modulationType, dataBits, ~, ~, total_sc, preamble, ~, window_size, preamble_cp_size, toa_error_estimate, showPlot)
%% OFDM Frame Synchronization, Channel Estimation & BER Simulation Kernel
% Course: ELEC-H422 / H423 (Wireless Communications & Propagation)
%
% Passes a time-domain OFDM burst (preamble + data blocks) through multipath CIR h
% and AWGN, performs Schmidl-Cox-style sliding autocorrelation ToA estimation,
% FFT demodulation, zero-forcing preamble channel estimation, and computes both
% unequalized and equalized BER along with CTF MSE.

if nargin < 21 || isempty(showPlot)
    showPlot = false;
end

%% Convolution with Multipath Channel & AWGN
timeSignalRx = conv(timeSignalTx, h);
power_sig_time = mean(abs(timeSignalRx).^2);
power_noise_time = power_sig_time / (10^(snr/10));
Noise_signal_time = sqrt(power_noise_time / 2) * (randn(1, size(timeSignalRx, 2)) + 1j * randn(1, size(timeSignalRx, 2)));
timeSignalRx_noisy = timeSignalRx + Noise_signal_time;

%% Time-Domain Synchronization (Sliding Autocorrelation)
num_samples = length(timeSignalRx_noisy);
num_windows = num_samples - window_size;
A = sum(timeSignalRx_noisy(1:total_sc) .* conj(timeSignalRx_noisy(1 + total_sc:2*total_sc)));
E = sum(abs(timeSignalRx_noisy(1:2*total_sc)).^2);

auto_corr = zeros(1, num_windows);
for n = 1:num_windows
    A = A - timeSignalRx_noisy(n) * conj(timeSignalRx_noisy(n + total_sc)) + timeSignalRx_noisy(n + total_sc) * conj(timeSignalRx_noisy(n + window_size));
    E = E - abs(timeSignalRx_noisy(n))^2 + abs(timeSignalRx_noisy(n + window_size))^2;
    auto_corr(n) = abs(A) / (E + eps);
end

if snr <= 5
    auto_corr = movmean(auto_corr, 15);
end

autocorr_values_avg_smooth = movmean(abs(auto_corr), preamble_cp_size);
[~, toa_estimate_smooth] = max(autocorr_values_avg_smooth);
toa_estimate_smooth = toa_estimate_smooth - size_cp;

ToA_error = (toa_estimate_smooth - true_ToA)^2;

% Shift back slightly inside cyclic prefix to absorb negative timing jitter
toa_index_corr = toa_estimate_smooth - toa_error_estimate;
if toa_index_corr + size_Signaltx > length(timeSignalRx_noisy)
    timeSignalRxCorr = timeSignalRx_noisy(1, toa_index_corr:end);
elseif toa_index_corr <= 0
    timeSignalRxCorr = timeSignalRx_noisy(1, 1:size_Signaltx-1);
else
    timeSignalRxCorr = timeSignalRx_noisy(1, toa_index_corr:toa_index_corr + size_Signaltx - 1);
end

num_samples = length(timeSignalRxCorr);
mod_value = mod(num_samples, size_cp + total_sc);
if mod_value ~= 0
    timeSignalRxCorr = timeSignalRxCorr(1:end - mod_value);
end

if showPlot
    figure;
    plot(abs(auto_corr)); hold on;
    plot(autocorr_values_avg_smooth, 'LineWidth', 1.4); hold off;
    title(sprintf('Schmidl-Cox Autocorrelation Timing Plateau (SNR = %d dB)', snr));
    xlabel('Sample Index'); ylabel('Normalized Metric');
    legend('Raw Metric', 'Smoothed (2\times CP)', 'Location', 'best'); grid on;
end

%% CP Removal & FFT Demodulation
preamble_rx_t = timeSignalRxCorr(1:2*(total_sc + size_cp));
preamble_rx_t(:, 1:preamble_cp_size) = [];
preamble_rx_t = reshape(preamble_rx_t, [], 2).';

timeBlocksRx = reshape(timeSignalRxCorr(1, 2*(total_sc + size_cp) + 1:end), total_sc + size_cp, []).';
timeBlocksRx(:, 1:size_cp, :) = [];
timeBlocksRx = [preamble_rx_t; timeBlocksRx];
OFDMBlocksRx = fft(timeBlocksRx, [], 2);

%% Preamble-Based Channel Estimation & CTF MSE
H_est = OFDMBlocksRx(1, :) ./ preamble;
H_est(isnan(H_est) | isinf(H_est)) = 0;
mse = mean(abs(abs(H(H_est ~= 0)) - abs(H_est(H_est ~= 0))).^2);

%% Zero-Forcing Frequency-Domain Equalization & BER
OFDMBlocksTxEstimate = OFDMBlocksRx(3:end, :) ./ H_est;
OFDMBlocksTxEstimate(isnan(OFDMBlocksTxEstimate) | isinf(OFDMBlocksTxEstimate)) = 0;

dataRxNotEqualized = reshape((OFDMBlocksRx(3:end, dataIndices + num_0_sc)).', [], 1);
demappedDataNotEqualized = demapping(dataRxNotEqualized, data_nbps, modulationType(1), 1);
error_neq = norm(dataBits(1:length(demappedDataNotEqualized)) - demappedDataNotEqualized, 1);
BER_notEqualized = error_neq / length(dataBits);

dataRxEqualized = reshape((OFDMBlocksTxEstimate(:, dataIndices + num_0_sc)).', [], 1);
demappedDataEqualized = demapping(dataRxEqualized, data_nbps, modulationType(1), 1);
error_eq = norm(dataBits(1:length(demappedDataEqualized)) - demappedDataEqualized, 1);
BER_Equalized = error_eq / length(dataBits);
end
