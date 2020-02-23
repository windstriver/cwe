U_ref = 11.1438;
Iu_ref = 0.116;
% sigma_u_ref = U_ref * Iu_ref;
sigma_u_ref = 1;
Lux = 0.5;

Suu_1d = @(omega) sigma_u_ref^2 * 2*Lux/(2*pi*U_ref) .* (1+70.8*(omega*Lux/(2*pi*U_ref)).^2).^(-5/6);

% integral(@(omega) Suu_1d(omega), -Inf, Inf)

% Coherency decay coefficient
Cuy = 10;
Cuz = 10;
% Coh_u = @(omega, y, z) exp(-Cu*(abs(y)+abs(z)).*abs(omega)./(2*pi*U_ref));

% Spectra of u
% Suu = @(k1, k2, k3) sigma_u_ref^2 * 2*Lux/(2*pi*U_ref) .* (1+70.8*(k1*Lux/(2*pi*U_ref)).^2).^(-5/6) .* ...
%     (1/pi)*(Cuy*abs(k1)/(2*pi*U_ref))./((Cuy*abs(k1)/(2*pi*U_ref)).^2+k2.^2) .* ...
%     (1/pi)*(Cuz*abs(k1)/(2*pi*U_ref))./((Cuz*abs(k1)/(2*pi*U_ref)).^2+k3.^2);

Suu = @(k1, k2) sigma_u_ref^2 * 2*Lux/(2*pi*U_ref) .* (1+70.8*(k1*Lux/(2*pi*U_ref)).^2).^(-5/6) .* ...
    (1./pi).*(Cuy.*abs(k1)./(2*pi*U_ref)) ./ ((Cuy.*abs(k1)./(2*pi*U_ref)).^2+k2.^2);

%% inlet plane
% Cut-off wavenumbers
k1u = 1000 * (2*pi);
k2u = 500 * (2*pi);
epsilon = 1 - integral2(Suu, 0, k1u, 0, k2u) / integral2(Suu, 0, Inf, 0, Inf);
fprintf('Relative truncation error with k1u=%6i and k2u=%6i is %5.3f. \n', k1u, k2u, epsilon);

% Discretization of wavenumbers
N1 = 10000;
N2 = 8000;
delta_k1 = k1u / N1;
delta_k2 = k2u / N2;
Lx10 = 2*pi/delta_k1;
Lx20 = 2*pi/delta_k2;

% Generation parameter
M1 = 2^(round(log(N1)/log(2))+2);
M2 = 2^(round(log(N2)/log(2))+2);
% M3 = 2^8;
delta_x1 = 2*pi/(M1*delta_k1);
delta_x2 = 2*pi/(M2*delta_k2);

k1_vec = (0:(M1-1))*delta_k1;
k2_vec = (0:(M2-1))*delta_k2;
% k3_vec = (0:(M3-1))*delta_k3;

% [k1_grid, k2_grid, k3_grid] = ndgrid(k1_vec, k2_vec, k3_vec);
[k1_grid, k2_grid] = ndgrid(k1_vec, k2_vec);

% A_grid = sqrt(2*arrayfun(Suu, k1_grid, k2_grid)*delta_k1*delta_k2);
A_grid = sqrt(2*Suu(k1_grid, k2_grid)*delta_k1*delta_k2);

% A_grid(:,:,1) = 0;
% A_grid(:,1,:) = 0;
% A_grid(1,:,:) = 0;
A_grid(:,1) = 0; A_grid(1,:) = 0;

2*sum(sum(A_grid.^2))

%% Generation
% Phi1 = 2*pi*rand(M1, M2);
% Phi2 = 2*pi*rand(M1, M2);

% B1 = sqrt(2) * A_grid .* exp(1i*Phi1);
% B2 = sqrt(2) * A_grid .* exp(1i*Phi2);

% u = real((M1*M2) * ifft2(B1) + M1*ifft(fft(B2, [], 2), [], 1));
