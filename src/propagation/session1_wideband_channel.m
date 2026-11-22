%% Wireless Propagation & MIMO-OFDM — Session 1: Wideband Street-Canyon Channel
% Course: ELEC-H422 / H423 (Wireless Communications & Propagation)
%
% Models wideband LOS/NLOS channel transfer functions (CTF), frequency correlation,
% coherence bandwidth, channel impulse response (CIR), power delay profile (PDP),
% and RMS delay spread across B in [1, 10, 100] MHz and reflection orders rho in [0, 2, 10].
% Requires proprietary urban-canyon models in ../../functions/*.p.

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(rootDir, 'functions'));
addpath(fullfile(rootDir, 'src', 'common'));

%% Channel Parameters
ch.Ptx = 0.1;
ch.wt = 8;
ch.wb = -8;
ch.rho = 5;
ch.fc = 3.6e9;
ch.B = 100e6;
ch.Q = 1024;
c = 3e8;
lambda_fc = c / ch.fc;
dr = lambda_fc / 4;
ch.beta = 2*pi / lambda_fc;
ch.A = ch.Ptx * lambda_fc^2 / (4*pi)^2;
ch.path_params.Npaths = 5;

%% Wideband Channel Transfer Function Sweep
B_list = [1e6 10e6 100e6];
rho_list = [0 2 10];
Rx_pos_x = 12:0.1:13;
Rx_pos_y = 6:0.1:7;
H_los_stack = zeros(length(rho_list), length(B_list), length(Rx_pos_y), ch.Q);
H_nlos_stack = zeros(length(rho_list), length(B_list), length(Rx_pos_y), ch.Q);

fig_los = figure; fig_los.Tag = 'LosCTF';
fig_nlos = figure; fig_nlos.Tag = 'nLosCTF';

for nb_reflex = 1:length(rho_list)
    ch.rho = rho_list(nb_reflex);
    for i = 1:length(B_list)
        ch.B = B_list(i);
        delta_f = ch.B / ch.Q;
        SCIndices = -ch.Q/2 + 1 : ch.Q/2;
        freqSC = SCIndices * delta_f;
        for j = 1:length(Rx_pos_x)
            ch.Rx_pos_x = Rx_pos_x(j);
            ch.Rx_pos_y = Rx_pos_y(j);
            [H_los, H_nlos, ~] = getWideBand_open(ch);
            H_los_stack(nb_reflex, i, j, :) = H_los;
            H_nlos_stack(nb_reflex, i, j, :) = H_nlos;
        end

        figure(fig_los);
        subplot(3, 3, (nb_reflex-1)*3 + i);
        plot(freqSC, 10*log10(abs(squeeze(H_los_stack(nb_reflex, i, :, :)))).');
        title(sprintf('CTF LOS %dMHz BW, %d Reflex', B_list(i)/1e6, ch.rho));
        xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]'); grid on;

        if nb_reflex > 1
            figure(fig_nlos);
            subplot(2, 3, (nb_reflex-2)*3 + i);
            plot(freqSC, 10*log10(abs(squeeze(H_nlos_stack(nb_reflex, i, :, :)))).');
            title(sprintf('CTF NLOS %dMHz BW, %d Reflex', B_list(i)/1e6, ch.rho));
            xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]'); grid on;
        end
    end
end

%% Frequency Correlation & Coherence Bandwidth
H_select_LOS = squeeze(H_los_stack(:, :, 10, :));
H_select_NLOS = squeeze(H_nlos_stack(:, :, 10, :));

nb_real = 250;
corr_los = zeros(length(rho_list), length(B_list), nb_real + 1);
corr_nlos = zeros(length(rho_list), length(B_list), nb_real + 1);
for i = 0:nb_real
    delta_f = i;
    corr_los(:, :, i+1) = mean(H_select_LOS(:, :, delta_f+1:end) .* conj(H_select_LOS(:, :, 1:end-delta_f)), 3);
    corr_nlos(:, :, i+1) = mean(H_select_NLOS(:, :, delta_f+1:end) .* conj(H_select_NLOS(:, :, 1:end-delta_f)), 3);
end

corr_normalized_los = corr_los ./ corr_los(:, :, 1);
corr_normalized_nlos = corr_nlos ./ corr_nlos(:, :, 1);

threshold = 0.7;
coherence_bw_los = zeros(3, 3);
coherence_bw_nlos = zeros(3, 3);
fig_FreqCorr = figure; fig_FreqCorr.Tag = 'FreqCorr';

for i = 1:size(corr_normalized_los, 1)
    ch.rho = rho_list(i);
    for j = 1:size(corr_normalized_los, 2)
        ch.B = B_list(j);
        delta_f = ch.B / ch.Q;
        SCIndices = -ch.Q/2 + 1 : ch.Q/2;
        freqSC = ch.fc + SCIndices * delta_f;

        idx_los = find(abs(corr_normalized_los(i, j, :)) <= threshold, 1);
        idx_nlos = find(abs(corr_normalized_nlos(i, j, :)) <= threshold, 1);
        if ~isempty(idx_los), coherence_bw_los(i, j) = freqSC(idx_los) - freqSC(1); else, coherence_bw_los(i, j) = ch.B; end
        if ~isempty(idx_nlos), coherence_bw_nlos(i, j) = freqSC(idx_nlos) - freqSC(1); else, coherence_bw_nlos(i, j) = ch.B; end

        figure(fig_FreqCorr);
        subplot(3, 3, (i-1)*3 + j);
        plot(abs(squeeze(corr_normalized_los(i, j, :))).'); hold on;
        if i > 1
            plot(abs(squeeze(corr_normalized_nlos(i, j, :))).');
            legend('LOS', 'NLOS', 'Location', 'best');
        end
        title(sprintf('Freq Corr %dMHz, %d Reflex', B_list(j)/1e6, ch.rho));
        xlabel('\Delta f bin'); ylabel('|R_H|'); grid on;
    end
end

%% Channel Impulse Response, PDP & Delay Spread
h_t_los = ifft(H_select_LOS, [], 3);
h_t_nlos = ifft(H_select_NLOS, [], 3);
Ts = (1 ./ B_list).';
size_crop = 100;
tau = repmat(permute((0:size_crop-1) .* Ts, [3 1 2]), 3, 1, 1);

h_t_los = h_t_los(:, :, 1:size_crop) ./ max(abs(h_t_los(:, :, 1:size_crop)), [], 3);
h_t_nlos = h_t_nlos(:, :, 1:size_crop) ./ max(abs(h_t_nlos(:, :, 1:size_crop)), [], 3);
PDP_los = abs(h_t_los).^2;
PDP_nlos = abs(h_t_nlos).^2;
PDP_los = PDP_los ./ max(PDP_los, [], 3);
PDP_nlos = PDP_nlos ./ max(PDP_nlos, [], 3);

P_los = sum(PDP_los, 3);
P_nlos = sum(PDP_nlos, 3);
tau_mean_los = (1 ./ P_los) .* sum(PDP_los .* tau, 3);
tau_mean_nlos = (1 ./ P_nlos) .* sum(PDP_nlos .* tau, 3);
sigma_tau_los = sqrt((1 ./ P_los) .* sum(PDP_los .* tau.^2, 3) - tau_mean_los.^2);
sigma_tau_nlos = sqrt((1 ./ P_nlos) .* sum(PDP_nlos .* tau.^2, 3) - tau_mean_nlos.^2);

fig_h_los = figure; fig_h_nlos = figure;
fig_PDP_los = figure; fig_PDP_nlos = figure;
for i = 1:3
    for j = 1:3
        figure(fig_h_los); subplot(3, 3, (i-1)*3 + j);
        stem(abs(squeeze(h_t_los(i, j, :)))); title(sprintf('CIR LOS %dMHz, %d Reflex', B_list(j)/1e6, rho_list(i)));
        figure(fig_PDP_los); subplot(3, 3, (i-1)*3 + j);
        stem(abs(squeeze(PDP_los(i, j, :)))); title(sprintf('PDP LOS %dMHz, %d Reflex', B_list(j)/1e6, rho_list(i)));
        figure(fig_h_nlos); subplot(3, 3, (i-1)*3 + j);
        stem(abs(squeeze(h_t_nlos(i, j, :)))); title(sprintf('CIR NLOS %dMHz, %d Reflex', B_list(j)/1e6, rho_list(i)));
        figure(fig_PDP_nlos); subplot(3, 3, (i-1)*3 + j);
        stem(abs(squeeze(PDP_nlos(i, j, :)))); title(sprintf('PDP NLOS %dMHz, %d Reflex', B_list(j)/1e6, rho_list(i)));
    end
end

fprintf('Session 1 Summary | Coherence BW LOS (MHz, rho=10): %s\n', sprintf('%.2f ', coherence_bw_los(3,:)/1e6));
fprintf('Session 1 Summary | Delay spread sigma_tau LOS (ns, rho=10): %s\n', sprintf('%.1f ', sigma_tau_los(3,:)*1e9));

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
        exportgraphics(figs(fi), fullfile(pngDir, sprintf('session1_wb_fig%02d.png', fi)), 'Resolution', 450, 'BackgroundColor', 'white');
        exportgraphics(figs(fi), fullfile(pdfDir, sprintf('session1_wb_fig%02d.pdf', fi)), 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
