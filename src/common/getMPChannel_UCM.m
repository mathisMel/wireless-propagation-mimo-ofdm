function [H, stochastic_RandG] = getMPChannel_UCM(path_params, Ptx, lambda, B, Q)

Rshift = path_params.shift_paths;
Npaths = path_params.Npaths;
Rmax_init = path_params.init_max;

% Rshift = 10; % shift MPC distances
% Rmax_init = 80; % initial max distance
D = sort((randperm(ceil(Rmax_init)*10,Npaths)./10)+Rshift, 'ascend');
D = D + randn(size(D));
dq = B/Q;
q = 0:Q-1;

A = sqrt(Ptx*lambda^2./(4*pi).^2);

Gamma = getGamma_UCM(randi([30, 45], size(D)), 3, ones(size(D)));

% Frequency-Delay Duality: Deterministic approachs
stochastic_RandG = [D; Gamma].';
tau_p = D ./ 3e8;
H = sum(A .* Gamma .* exp(1j*randi([1 314], 1, Npaths)./100) .* exp(-1j*2*pi*q.'.*tau_p*dq)./D , 2);

end

function Gamma = getGamma_UCM(Phi, eps_r, reflection_order)
Gamma_numerator = (cosd(Phi) - sqrt(eps_r).*sqrt(1-(1/eps_r)*sind(Phi).^2));
Gamma_denominator = (cosd(Phi) + sqrt(eps_r).*sqrt(1-(1/eps_r)*sind(Phi).^2));

Gamma_vertical_intermediate = abs(Gamma_numerator ./ Gamma_denominator); % abs or not?


Gamma = Gamma_vertical_intermediate.^reflection_order;
end