function [H_los, H_nlos, RandTheta] = getWideBand_open(ch)

Rx_pos_x = ch.Rx_pos_x;
Rx_pos_y = ch.Rx_pos_y;
A = ch.A;
wt = ch.wt;
wb = abs(ch.wb);
rho = ch.rho;
Ptx = ch.Ptx;

B = ch.B;
T = 1/B;
Q = ch.Q;
fc = ch.fc;
lambda_fc = 3e8/fc;
dq = 1/(Q*T);
subQ = (0:Q-1)*dq;
lambda_q = 3e8./subQ;
% lambda_q = lambda_q + lambda_fc;
beta_q = 2*pi./lambda_q;
%%

Rx = [Rx_pos_x Rx_pos_y];

rayInfo = ucm_getRayInfo(Rx, [wt wb], rho, 0);
% H_los = zeros(Q, 1);
% H_nlos = zeros(Q, 1);
ch.path_params.shift_paths = round(max(mink(rayInfo.R, 2)));

RandG_UCM = sortrows([rayInfo.R rayInfo.Gamma]);

if rho ~= 0
    RandG = zeros(ch.path_params.Npaths, size(RandG_UCM, 1)-1, 2);
    for r = 1:size(RandG_UCM, 1)-1

        ch.path_params.shift_paths = RandG_UCM(r, 1);
        ch.path_params.init_max = RandG_UCM(r+1, 1);
        [~, RandG(:, r, :)]= getMPChannel_UCM(ch.path_params, Ptx, lambda_fc, B, Q);
        RandG(1, r, :) = RandG_UCM(r, :);
    end
    RandG = reshape(RandG, ch.path_params.Npaths*(size(RandG_UCM, 1)-1), 2);
    % Rmax = ch.path_params.max_delay*3e8;
    % RandG((RandG(:, 1) > Rmax), :) = [];
else
    RandG = RandG_UCM;
end

q_vec = 0:Q-1;
random_phases = randi([1 360], [size(RandG, 1), 1])./180*pi;
H_temp = exp(-1j*2*pi*q_vec.*RandG(:,1)./(3e8*Q*T)).*RandG(:,2)./RandG(:,1) .* exp(1j*random_phases);
H_los = sqrt(A)*sum(H_temp, 1);
H_nlos = sqrt(A)*sum(H_temp(2:end, :), 1);

RandTheta = sortrows([rayInfo.R rayInfo.AoD]);

end