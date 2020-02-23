%% Parameters
h5temp = 'temp.h5';
h5file = 'inflowTurb.h5';

% wind velocity statistics
U_ref = 11.1438;
Iu_ref = 0.116;
% sigma_u_ref = U_ref * Iu_ref;
sigma_u_ref = 1;
sigma_v_ref = 0.78*sigma_u_ref;
sigma_w_ref = 0.55*sigma_u_ref;
Lux = 0.5;
Lvx = 0.5 * Lux;
Lwx = 0.5 * Lux;

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

%% Spectral matrix
% Cut-off wavenumbers
k1u = 400 * (2*pi);
k2u = 40 * (2*pi);
k3u = 40 * (2*pi);
% epsilon = 1 - sqrt(integral3(Suu, 0, k1u, 0, k2u, 0, k3u) / integral3(Suu, 0, Inf, 0, Inf, 0, Inf));
% fprintf('Relative truncation error with k1u=%6i and k2u=%6i is %5.3f. \n', k1u, k2u, epsilon);

% Discretization of wavenumbers
N1 = 4800;
N2 = 400;
N3 = 400;
delta_k1 = k1u / N1;
delta_k2 = k2u / N2;
delta_k3 = k3u / N3;
Lx10 = 2*pi/delta_k1;
Lx20 = 2*pi/delta_k2;
Lx30 = 2*pi/delta_k3;

M1 = 2^(round(log(N1)/log(2))+2);
M2 = 2^(round(log(N2)/log(2))+2);
M3 = 2^(round(log(N3)/log(2))+2);

delta_x1 = 2*pi/(M1*delta_k1);
delta_x2 = 2*pi/(M2*delta_k2);
delta_x3 = 2*pi/(M3*delta_k3);

k1_vec = (0:(M1-1))*delta_k1;
k2_vec = (0:(M2-1))*delta_k2;
k3_vec = (0:(M3-1))*delta_k3;

% Spectral matrix
try
    h5create(h5temp, '/Au', [M1, M2, M3]);
    h5create(h5temp, '/Av', [M1, M2, M3]);
    h5create(h5temp, '/Aw', [M1, M2, M3]);
catch
end

disp('Generate spectral matrix:');
for k = 1:M3
    [k1_grid, k2_grid, k3_grid] = ndgrid(k1_vec, k2_vec, k3_vec(k));
    Au = sqrt(2*Suu(k1_grid, k2_grid, k3_grid)*delta_k1*delta_k2*delta_k3);
    Av = sqrt(2*Svv(k1_grid, k2_grid, k3_grid)*delta_k1*delta_k2*delta_k3);
    Aw = sqrt(2*Sww(k1_grid, k2_grid, k3_grid)*delta_k1*delta_k2*delta_k3);

    h5write(h5temp, '/Au', Au, [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Av', Av, [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Aw', Aw, [1 1 k], [M1 M2 1]);
end

h5write(h5temp, '/Au', zeros(M1,M2), [1 1 1], [M1 M2 1]);
h5write(h5temp, '/Au', zeros(M1,1,M3), [1 1 1], [M1 1 M3]);
h5write(h5temp, '/Au', zeros(1,M2,M3), [1 1 1], [1 M2 M3]);

h5write(h5temp, '/Av', zeros(M1,M2), [1 1 1], [M1 M2 1]);
h5write(h5temp, '/Av', zeros(M1,1,M3), [1 1 1], [M1 1 M3]);
h5write(h5temp, '/Av', zeros(1,M2,M3), [1 1 1], [1 M2 M3]);

h5write(h5temp, '/Aw', zeros(M1,M2), [1 1 1], [M1 M2 1]);
h5write(h5temp, '/Aw', zeros(M1,1,M3), [1 1 1], [M1 1 M3]);
h5write(h5temp, '/Aw', zeros(1,M2,M3), [1 1 1], [1 M2 M3]);

disp('Finished gnerating spectral matrix');

%% Generation
% Create HDF5 dataset
try
    h5create(h5temp, '/Phi1', [M1 M2 M3]);
    h5create(h5temp, '/Phi2', [M1 M2 M3]);
    h5create(h5temp, '/Phi3', [M1 M2 M3]);
    h5create(h5temp, '/Phi4', [M1 M2 M3]);
    h5create(h5temp, '/re1', [M1 M2 M3]);
    h5create(h5temp, '/im1', [M1 M2 M3]);
    h5create(h5temp, '/re2', [M1 M2 M3]);
    h5create(h5temp, '/im2', [M1 M2 M3]);
    h5create(h5temp, '/re3', [M1 M2 M3]);
    h5create(h5temp, '/im3', [M1 M2 M3]);
    h5create(h5temp, '/re4', [M1 M2 M3]);
    h5create(h5temp, '/im4', [M1 M2 M3]);
catch
end

try
    h5create(h5file, '/gamma_uw', [M2 M3]);
    h5create(h5file, '/u', [M1 M2 M3]);
    h5create(h5file, '/w1', [M1 M2 M3]);
    h5create(h5file, '/v', [M1 M2 M3]);
    h5create(h5file, '/w3', [M1 M2 M3]);
catch
end

%% First sample for u
disp('Start generating first sample for u:');
for k = 1:M3
    h5write(h5temp, '/Phi1', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi2', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi3', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi4', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
end

for k = 1:M3
    Au = h5read(h5temp, '/Au', [1 1 k], [M1 M2 1]);

    Phi1 = h5read(h5temp, '/Phi1', [1 1 k], [M1 M2 1]);
    Phi2 = h5read(h5temp, '/Phi2', [1 1 k], [M1 M2 1]);
    Phi3 = h5read(h5temp, '/Phi3', [1 1 k], [M1 M2 1]);
    Phi4 = h5read(h5temp, '/Phi4', [1 1 k], [M1 M2 1]);
    
    u1 = (M1*M2) * ifft(ifft(sqrt(2)*Au.*exp(1i*Phi1), [], 2), [], 1);
    u2 = (M1*M2) * ifft(ifft(sqrt(2)*Au.*exp(1i*Phi2), [], 2), [], 1);
    u3 = (M1   ) * ifft( fft(sqrt(2)*Au.*exp(1i*Phi3), [], 2), [], 1);
    u4 = (M1   ) * ifft( fft(sqrt(2)*Au.*exp(1i*Phi4), [], 2), [], 1);
    
    h5write(h5temp, '/re1', real(u1), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im1', imag(u1), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re2', real(u2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im2', imag(u2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re3', real(u3), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im3', imag(u3), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re4', real(u4), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im4', imag(u4), [1 1 k], [M1 M2 1]);
end

for j = 1:M2
    re1 = h5read(h5temp, '/re1', [1 j 1], [M1 1 M3]);
    im1 = h5read(h5temp, '/im1', [1 j 1], [M1 1 M3]);
    re2 = h5read(h5temp, '/re2', [1 j 1], [M1 1 M3]);
    im2 = h5read(h5temp, '/im2', [1 j 1], [M1 1 M3]);
    re3 = h5read(h5temp, '/re3', [1 j 1], [M1 1 M3]);
    im3 = h5read(h5temp, '/im3', [1 j 1], [M1 1 M3]);
    re4 = h5read(h5temp, '/re4', [1 j 1], [M1 1 M3]);
    im4 = h5read(h5temp, '/im4', [1 j 1], [M1 1 M3]);
    
    u = M3 * ifft(re1+1i*im1, [], 3) + ...
              fft(re2+1i*im2, [], 3) + ...
        M3 * ifft(re3+1i*im3, [], 3) + ...
              fft(re4+1i*im4, [], 3);

    h5write(h5file, '/u', real(u), [1 j 1], [M1 1 M3]);
end

disp('Finish generating first sample for u.');

%% First sample for w
disp('Start generating first sample for w:');

for k = 1:M3
    Au = h5read(h5temp, '/Aw', [1 1 k], [M1 M2 1]);

    Phi1 = h5read(h5temp, '/Phi1', [1 1 k], [M1 M2 1]);
    Phi2 = h5read(h5temp, '/Phi2', [1 1 k], [M1 M2 1]);
    Phi3 = h5read(h5temp, '/Phi3', [1 1 k], [M1 M2 1]);
    Phi4 = h5read(h5temp, '/Phi4', [1 1 k], [M1 M2 1]);
    
    u1 = (M1*M2) * ifft(ifft(sqrt(2)*Au.*exp(1i*Phi1), [], 2), [], 1);
    u2 = (M1*M2) * ifft(ifft(sqrt(2)*Au.*exp(1i*Phi2), [], 2), [], 1);
    u3 = (M1   ) * ifft( fft(sqrt(2)*Au.*exp(1i*Phi3), [], 2), [], 1);
    u4 = (M1   ) * ifft( fft(sqrt(2)*Au.*exp(1i*Phi4), [], 2), [], 1);
    
    h5write(h5temp, '/re1', real(u1), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im1', imag(u1), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re2', real(u2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im2', imag(u2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re3', real(u3), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im3', imag(u3), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re4', real(u4), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im4', imag(u4), [1 1 k], [M1 M2 1]);
end

for j = 1:M2
    re1 = h5read(h5temp, '/re1', [1 j 1], [M1 1 M3]);
    im1 = h5read(h5temp, '/im1', [1 j 1], [M1 1 M3]);
    re2 = h5read(h5temp, '/re2', [1 j 1], [M1 1 M3]);
    im2 = h5read(h5temp, '/im2', [1 j 1], [M1 1 M3]);
    re3 = h5read(h5temp, '/re3', [1 j 1], [M1 1 M3]);
    im3 = h5read(h5temp, '/im3', [1 j 1], [M1 1 M3]);
    re4 = h5read(h5temp, '/re4', [1 j 1], [M1 1 M3]);
    im4 = h5read(h5temp, '/im4', [1 j 1], [M1 1 M3]);
    
    u = M3 * ifft(re1+1i*im1, [], 3) + ...
              fft(re2+1i*im2, [], 3) + ...
        M3 * ifft(re3+1i*im3, [], 3) + ...
              fft(re4+1i*im4, [], 3);

    h5write(h5file, '/w1', real(u), [1 j 1], [M1 1 M3]);
end

disp('Finish generating first sample for w.');

gamma_uw = zeros(M2, M3);
for j=1:M2
    for k = 1:M3
        u = h5read(h5file, '/u', [1 j k], [M1 1 1]);
        w = h5read(h5file, '/w1', [1 j k], [M1 1 1]);
        corrmat = corrcoef(u, w1);
        gamma_uw(j,k) = corrmat(1,2);
    end
end

h5write(h5file, '/gamma_uw', gamma_uw);

%% Second sample for v
disp('Start generating second sample for v:');

for k = 1:M3
    h5write(h5temp, '/Phi1', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi2', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi3', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi4', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
end

for k = 1:M3
    Au = h5read(h5temp, '/Av', [1 1 k], [M1 M2 1]);

    Phi1 = h5read(h5temp, '/Phi1', [1 1 k], [M1 M2 1]);
    Phi2 = h5read(h5temp, '/Phi2', [1 1 k], [M1 M2 1]);
    Phi3 = h5read(h5temp, '/Phi3', [1 1 k], [M1 M2 1]);
    Phi4 = h5read(h5temp, '/Phi4', [1 1 k], [M1 M2 1]);
    
    u1 = (M1*M2) * ifft(ifft(sqrt(2)*Au.*exp(1i*Phi1), [], 2), [], 1);
    u2 = (M1*M2) * ifft(ifft(sqrt(2)*Au.*exp(1i*Phi2), [], 2), [], 1);
    u3 = (M1   ) * ifft( fft(sqrt(2)*Au.*exp(1i*Phi3), [], 2), [], 1);
    u4 = (M1   ) * ifft( fft(sqrt(2)*Au.*exp(1i*Phi4), [], 2), [], 1);
    
    h5write(h5temp, '/re1', real(u1), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im1', imag(u1), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re2', real(u2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im2', imag(u2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re3', real(u3), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im3', imag(u3), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re4', real(u4), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im4', imag(u4), [1 1 k], [M1 M2 1]);
end

for j = 1:M2
    re1 = h5read(h5temp, '/re1', [1 j 1], [M1 1 M3]);
    im1 = h5read(h5temp, '/im1', [1 j 1], [M1 1 M3]);
    re2 = h5read(h5temp, '/re2', [1 j 1], [M1 1 M3]);
    im2 = h5read(h5temp, '/im2', [1 j 1], [M1 1 M3]);
    re3 = h5read(h5temp, '/re3', [1 j 1], [M1 1 M3]);
    im3 = h5read(h5temp, '/im3', [1 j 1], [M1 1 M3]);
    re4 = h5read(h5temp, '/re4', [1 j 1], [M1 1 M3]);
    im4 = h5read(h5temp, '/im4', [1 j 1], [M1 1 M3]);
    
    u = M3 * ifft(re1+1i*im1, [], 3) + ...
              fft(re2+1i*im2, [], 3) + ...
        M3 * ifft(re3+1i*im3, [], 3) + ...
              fft(re4+1i*im4, [], 3);

    h5write(h5file, '/v', real(u), [1 j 1], [M1 1 M3]);
end

disp('Finish generating second sample for v.');

%% Third sample for w
disp('Start generating Third sample for w:');

for k = 1:M3
    h5write(h5temp, '/Phi1', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi2', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi3', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/Phi4', 2*pi*rand(M1,M2), [1 1 k], [M1 M2 1]);
end

for k = 1:M3
    Au = h5read(h5temp, '/Aw', [1 1 k], [M1 M2 1]);

    Phi1 = h5read(h5temp, '/Phi1', [1 1 k], [M1 M2 1]);
    Phi2 = h5read(h5temp, '/Phi2', [1 1 k], [M1 M2 1]);
    Phi3 = h5read(h5temp, '/Phi3', [1 1 k], [M1 M2 1]);
    Phi4 = h5read(h5temp, '/Phi4', [1 1 k], [M1 M2 1]);
    
    u1 = (M1*M2) * ifft(ifft(sqrt(2)*Au.*exp(1i*Phi1), [], 2), [], 1);
    u2 = (M1*M2) * ifft(ifft(sqrt(2)*Au.*exp(1i*Phi2), [], 2), [], 1);
    u3 = (M1   ) * ifft( fft(sqrt(2)*Au.*exp(1i*Phi3), [], 2), [], 1);
    u4 = (M1   ) * ifft( fft(sqrt(2)*Au.*exp(1i*Phi4), [], 2), [], 1);
    
    h5write(h5temp, '/re1', real(u1), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im1', imag(u1), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re2', real(u2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im2', imag(u2), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re3', real(u3), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im3', imag(u3), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/re4', real(u4), [1 1 k], [M1 M2 1]);
    h5write(h5temp, '/im4', imag(u4), [1 1 k], [M1 M2 1]);
end

for j = 1:M2
    re1 = h5read(h5temp, '/re1', [1 j 1], [M1 1 M3]);
    im1 = h5read(h5temp, '/im1', [1 j 1], [M1 1 M3]);
    re2 = h5read(h5temp, '/re2', [1 j 1], [M1 1 M3]);
    im2 = h5read(h5temp, '/im2', [1 j 1], [M1 1 M3]);
    re3 = h5read(h5temp, '/re3', [1 j 1], [M1 1 M3]);
    im3 = h5read(h5temp, '/im3', [1 j 1], [M1 1 M3]);
    re4 = h5read(h5temp, '/re4', [1 j 1], [M1 1 M3]);
    im4 = h5read(h5temp, '/im4', [1 j 1], [M1 1 M3]);
    
    u = M3 * ifft(re1+1i*im1, [], 3) + ...
              fft(re2+1i*im2, [], 3) + ...
        M3 * ifft(re3+1i*im3, [], 3) + ...
              fft(re4+1i*im4, [], 3);

    h5write(h5file, '/w3', real(u), [1 j 1], [M1 1 M3]);
end

disp('Finish generating third sample for w.');

%% Save to HDF5
h5writeatt(h5file, '/', 'delta_x1', delta_x1);
h5writeatt(h5file, '/', 'M1', M1);
h5writeatt(h5file, '/', 'delta_x2', delta_x2);
h5writeatt(h5file, '/', 'M2', M2);
h5writeatt(h5file, '/', 'delta_x3', delta_x3);
h5writeatt(h5file, '/', 'M3', M3);
