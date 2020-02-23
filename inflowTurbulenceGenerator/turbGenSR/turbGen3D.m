U_ref = 11.1438;
Iu_ref = 0.116;
sigma_u_ref = U_ref * Iu_ref;
% sigma_u_ref = 1;
sigma_v_ref = 0.78*sigma_u_ref;
sigma_w_ref = 0.55*sigma_u_ref;
Lux = 0.5;
Lvx = 0.237 * Lux;
Lwx = 0.082 * Lux;

% Coherency decay coefficient
Cuy = 10;
Cuz = 10;

% Spectra
Suu = @(k1, k2, k3) sigma_u_ref^2 * 2*Lux/(2*pi*U_ref) .* ...
    (1+70.8*(k1*Lux/(2*pi*U_ref)).^2).^(-5/6) .* ...
    (1./pi).*(Cuy.*abs(k1)./(2*pi*U_ref)) ./ ((Cuy.*abs(k1)./(2*pi*U_ref)).^2+k2.^2) .* ...
    (1./pi).*(Cuz.*abs(k1)./(2*pi*U_ref)) ./ ((Cuz.*abs(k1)./(2*pi*U_ref)).^2+k3.^2);

Svv = @(k1, k2, k3) sigma_v_ref^2 * 2*Lvx/(2*pi*U_ref) .* ...
    (1+755.2*(k1*Lvx/(2*pi*U_ref)).^2) .* ...
    (1+283.2*(k1*Lvx/(2*pi*U_ref)).^2).^(-11/6) .* ...
    (1./pi).*(Cuy.*abs(k1)./(2*pi*U_ref)) ./ ((Cuy.*abs(k1)./(2*pi*U_ref)).^2+k2.^2) .* ...
    (1./pi).*(Cuz.*abs(k1)./(2*pi*U_ref)) ./ ((Cuz.*abs(k1)./(2*pi*U_ref)).^2+k3.^2);


Sww = @(k1, k2, k3) sigma_w_ref^2 * 2*Lwx/(2*pi*U_ref) .* ...
    (1+755.2*(k1*Lwx/(2*pi*U_ref)).^2) .* ...
    (1+283.2*(k1*Lwx/(2*pi*U_ref)).^2).^(-11/6) .* ...
    (1./pi).*(Cuy.*abs(k1)./(2*pi*U_ref)) ./ ((Cuy.*abs(k1)./(2*pi*U_ref)).^2+k2.^2) .* ...
    (1./pi).*(Cuz.*abs(k1)./(2*pi*U_ref)) ./ ((Cuz.*abs(k1)./(2*pi*U_ref)).^2+k3.^2);

%% inlet plane
% Cut-off wavenumbers
k1u = 400 * (2*pi);
k2u = 40 * (2*pi);
k3u = 40 * (2*pi);
epsilon = 1 - sqrt(integral3(Suu, 0, k1u, 0, k2u, 0, k3u) / integral3(Suu, 0, Inf, 0, Inf, 0, Inf));
fprintf('Relative truncation error with k1u=%6i and k2u=%6i is %5.3f. \n', k1u, k2u, epsilon);

% Discretization of wavenumbers
N1 = 4800;
N2 = 40*5;
N3 = 40*2;
delta_k1 = k1u / N1;
delta_k2 = k2u / N2;
delta_k3 = k3u / N3;
Lx10 = 2*pi/delta_k1;
Lx20 = 2*pi/delta_k2;
Lx30 = 2*pi/delta_k3;

% Generation parameter
M1 = 2^(round(log(N1)/log(2))+2);
M2 = 2^(round(log(N2)/log(2))+2);
M3 = 2^(round(log(N3)/log(2))+2);

delta_x1 = 2*pi/(M1*delta_k1);
delta_x2 = 2*pi/(M2*delta_k2);
delta_x3 = 2*pi/(M3*delta_k3);

k1_vec = (0:(M1-1))*delta_k1;
k2_vec = (0:(M2-1))*delta_k2;
k3_vec = (0:(M3-1))*delta_k3;

[k1_grid, k2_grid, k3_grid] = ndgrid(k1_vec, k2_vec, k3_vec);

Au_grid = sqrt(2*Suu(k1_grid, k2_grid, k3_grid)*delta_k1*delta_k2*delta_k3);
Av_grid = sqrt(2*Svv(k1_grid, k2_grid, k3_grid)*delta_k1*delta_k2*delta_k3);
Aw_grid = sqrt(2*Sww(k1_grid, k2_grid, k3_grid)*delta_k1*delta_k2*delta_k3);

Au_grid(:,:,1) = 0; Au_grid(:,1,:) = 0; Au_grid(1,:,:) = 0;
Av_grid(:,:,1) = 0; Av_grid(:,1,:) = 0; Av_grid(1,:,:) = 0;
Aw_grid(:,:,1) = 0; Aw_grid(:,1,:) = 0; Aw_grid(1,:,:) = 0;

sqrt(4*sum(sum(sum(Au_grid.^2))))

clear k1_vec k2_vec k3_vec
clear k1_grid k2_grid k3_grid

%% Generation
rho_uw = -0.2;

% First sample
Phi1 = 2*pi*rand(M1, M2, M3);
Phi2 = 2*pi*rand(M1, M2, M3);
Phi3 = 2*pi*rand(M1, M2, M3);
Phi4 = 2*pi*rand(M1, M2, M3);

u = real((M1*M2*M3) * ifft(ifft(ifft(sqrt(2) * Au_grid .* exp(1i*Phi1), [], 3), [], 2), [], 1) + ...
    (M1*M2) * ifft(ifft(fft(sqrt(2) * Au_grid .* exp(1i*Phi2), [], 3), [], 2), [], 1) + ...
    (M1*M3) * ifft(fft(ifft(sqrt(2) * Au_grid .* exp(1i*Phi3), [], 3), [], 2), [], 1) + ...
    (M1) * ifft(fft(fft(sqrt(2) * Au_grid .* exp(1i*Phi4), [], 3), [], 2), [], 1));

w1 = real((M1*M2*M3) * ifft(ifft(ifft(sqrt(2) * Aw_grid .* exp(1i*Phi1), [], 3), [], 2), [], 1) + ...
    (M1*M2) * ifft(ifft(fft(sqrt(2) * Aw_grid .* exp(1i*Phi2), [], 3), [], 2), [], 1) + ...
    (M1*M3) * ifft(fft(ifft(sqrt(2) * Aw_grid .* exp(1i*Phi3), [], 3), [], 2), [], 1) + ...
    (M1) * ifft(fft(fft(sqrt(2) * Aw_grid .* exp(1i*Phi4), [], 3), [], 2), [], 1));

gamma_uw = corrcoef(u, w1);
gamma_uw = gamma_uw(1,2);

% Second sample
Phi1 = 2*pi*rand(M1, M2, M3);
Phi2 = 2*pi*rand(M1, M2, M3);
Phi3 = 2*pi*rand(M1, M2, M3);
Phi4 = 2*pi*rand(M1, M2, M3);

v = real((M1*M2*M3) * ifft(ifft(ifft(sqrt(2) * Av_grid .* exp(1i*Phi1), [], 3), [], 2), [], 1) + ...
    (M1*M2) * ifft(ifft(fft(sqrt(2) * Av_grid .* exp(1i*Phi2), [], 3), [], 2), [], 1) + ...
    (M1*M3) * ifft(fft(ifft(sqrt(2) * Av_grid .* exp(1i*Phi3), [], 3), [], 2), [], 1) + ...
    (M1) * ifft(fft(fft(sqrt(2) * Av_grid .* exp(1i*Phi4), [], 3), [], 2), [], 1));

% Third sample
Phi1 = 2*pi*rand(M1, M2, M3);
Phi2 = 2*pi*rand(M1, M2, M3);
Phi3 = 2*pi*rand(M1, M2, M3);
Phi4 = 2*pi*rand(M1, M2, M3);

w3 = real((M1*M2*M3) * ifft(ifft(ifft(sqrt(2) * Aw_grid .* exp(1i*Phi1), [], 3), [], 2), [], 1) + ...
    (M1*M2) * ifft(ifft(fft(sqrt(2) * Aw_grid .* exp(1i*Phi2), [], 3), [], 2), [], 1) + ...
    (M1*M3) * ifft(fft(ifft(sqrt(2) * Aw_grid .* exp(1i*Phi3), [], 3), [], 2), [], 1) + ...
    (M1) * ifft(fft(fft(sqrt(2) * Aw_grid .* exp(1i*Phi4), [], 3), [], 2), [], 1));

w = rho_uw ./ gamma_uw .* w1 + sqrt(1-(rho_uw./gamma_uw).^2) .* w3;

%% Save to HDF5
h5file = 'inflowTurb.h5';

h5create(h5file, '/u', size(u));
h5write(h5file, '/u', u);
h5create(h5file, '/v', size(v));
h5write(h5file, '/v', v);
h5create(h5file, '/w', size(w));
h5write(h5file, '/w', w);

h5writeatt(h5file, '/', 'delta_x1', delta_x1);
h5writeatt(h5file, '/', 'M1', M1);
h5writeatt(h5file, '/', 'delta_x2', delta_x2);
h5writeatt(h5file, '/', 'M2', M2);
h5writeatt(h5file, '/', 'delta_x3', delta_x3);
h5writeatt(h5file, '/', 'M3', M3);
