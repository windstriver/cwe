clear; clc;
%% Spectral tensor parameter
% Kaimal model
gamma = 3.9;
z = 40; %[m]
L = 0.59 * z; %[m]
ae23 = 0.86; %[m^(4/3) s^(-1)

u_tau = sqrt(ae23 * z.^(2/3) / 3.2);

%% Simulation box parameters
Lx10 =  64 * L;
Lx20 =  16 * L;
Lx30 =  16 * L;

N1 = 2048;
N2 =  256;
N3 =  256;

delta_x1 = Lx10 / N1;
delta_x2 = Lx20 / N2;
delta_x3 = Lx30 / N3;

delta_k1 = 2*pi / Lx10;
delta_k2 = 2*pi / Lx20;
delta_k3 = 2*pi / Lx30;

%% Random phases
h5temp = 'turbGenTemp.h5';

if ~exist(h5temp, 'file')
    h5create(h5temp, '/Phiu', [N1 N2 N3], 'FillValue', 0);
    h5create(h5temp, '/Phiv', [N1 N2 N3], 'FillValue', 0);
    h5create(h5temp, '/Phiw', [N1 N2 N3], 'FillValue', 0);
end

for i = 2:N1/2
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiu', Phiu, [i 2 2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiu', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i N2/2+2 N3/2+2], [1 N2/2-1 N3/2-1]);
    
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiu', Phiu, [i N2/2+2 2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiu', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i 2 N3/2+2], [1 N2/2-1 N3/2-1]);

    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiu', Phiu, [i 2 N3/2+2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiu', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i N2/2+2 2], [1 N2/2-1 N3/2-1]);
    
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiu', Phiu, [i N2/2+2 N3/2+2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiu', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i 2 2], [1 N2/2-1 N3/2-1]);
    
    fprintf('Random phase for u generated at first index of %u\n', i);
end

for i = 2:N1/2
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiv', Phiu, [i 2 2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiv', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i N2/2+2 N3/2+2], [1 N2/2-1 N3/2-1]);
    
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiv', Phiu, [i N2/2+2 2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiv', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i 2 N3/2+2], [1 N2/2-1 N3/2-1]);

    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiv', Phiu, [i 2 N3/2+2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiv', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i N2/2+2 2], [1 N2/2-1 N3/2-1]);
    
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiv', Phiu, [i N2/2+2 N3/2+2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiv', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i 2 2], [1 N2/2-1 N3/2-1]);
    
    fprintf('Random phase for v generated at first index of %u\n', i);
end

for i = 2:N1/2
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiw', Phiu, [i 2 2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiw', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i N2/2+2 N3/2+2], [1 N2/2-1 N3/2-1]);
    
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiw', Phiu, [i N2/2+2 2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiw', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i 2 N3/2+2], [1 N2/2-1 N3/2-1]);

    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiw', Phiu, [i 2 N3/2+2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiw', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i N2/2+2 2], [1 N2/2-1 N3/2-1]);
    
    Phiu = 2*pi*rand(1, N2/2-1, N3/2-1);
    h5write(h5temp, '/Phiw', Phiu, [i N2/2+2 N3/2+2], [1 N2/2-1 N3/2-1]);
    h5write(h5temp, '/Phiw', -Phiu(1,end:-1:1,end:-1:1), [N1+2-i 2 2], [1 N2/2-1 N3/2-1]);
    
    fprintf('Random phase for w generated at first index of %u\n', i);
end

%% Cholesky decomposition of spectral tensor
h5C = 'C.h5';
if ~exist(h5C,'file')
    h5create(h5C, '/C11', [N1 N2 N3], 'FillValue', 0);
    h5create(h5C, '/C21', [N1 N2 N3], 'FillValue', 0);
    h5create(h5C, '/C31', [N1 N2 N3], 'FillValue', 0);
    h5create(h5C, '/C22', [N1 N2 N3], 'FillValue', 0);
    h5create(h5C, '/C32', [N1 N2 N3], 'FillValue', 0);
    h5create(h5C, '/C33', [N1 N2 N3], 'FillValue', 0);
end

k1_vec = [0, 1:(N1/2), -(N1/2-1):(-1)] * delta_k1;
k2_vec = [0, 1:(N2/2), -(N2/2-1):(-1)] * delta_k2;
k3_vec = [0, 1:(N3/2), -(N3/2-1):(-1)] * delta_k3;

for k = 1:N3
    [k1_grid, k2_grid, k3_grid] = ndgrid(k1_vec, k2_vec, k3_vec(k));

    % Using Eq. (47)
    % [RHS11, RHS22, RHS33, RHS12, RHS13, RHS23] = ...
    %     coefFS(k1_grid, k2_grid, k3_grid, Lx10, Lx20, Lx30, gamma, L, ae23, 21);

    % Using Eq. (46)
    [RHS11, RHS22, RHS33, RHS12, RHS13, RHS23] = ...
        spectralTensor(k1_grid, k2_grid, k3_grid, gamma, L, ae23);

    RHS11 = (2*pi)^3 / (Lx10*Lx20*Lx30) * RHS11;
    RHS22 = (2*pi)^3 / (Lx10*Lx20*Lx30) * RHS22;
    RHS33 = (2*pi)^3 / (Lx10*Lx20*Lx30) * RHS33;
    RHS12 = (2*pi)^3 / (Lx10*Lx20*Lx30) * RHS12;
    RHS13 = (2*pi)^3 / (Lx10*Lx20*Lx30) * RHS13;
    RHS23 = (2*pi)^3 / (Lx10*Lx20*Lx30) * RHS23;

    % Cholesky decomposition
    C11 = sqrt(RHS11);
    C21 = RHS12 ./ C11;
    C31 = RHS13 ./ C11;
    C22 = sqrt(RHS22 - C21.^2);
    C32 = (RHS23 - C21 .* C31) ./ C22;
    C33 = sqrt(RHS33 + 1e-12 - C32 .^2 - C31 .^2);

    if (k == 1) || (k == N2/2+1)
        C11 = zeros(N1, N2);
        C21 = zeros(N1, N2);
        C31 = zeros(N1, N2);
        C22 = zeros(N1, N2);
        C32 = zeros(N1, N2);
        C33 = zeros(N1, N2);
    end

    C11(1,:) = 0; C11(:,1) = 0;
    C11(N1/2+1,:) = 0; C11(:,N2/2+1) = 0;

    C21(1,:) = 0; C21(:,1) = 0;
    C21(N1/2+1,:) = 0; C21(:,N2/2+1) = 0;
    
    C31(1,:) = 0; C31(:,1) = 0;
    C31(N1/2+1,:) = 0; C31(:,N2/2+1) = 0;
    
    C22(1,:) = 0; C22(:,1) = 0;
    C22(N1/2+1,:) = 0; C22(:,N2/2+1) = 0;
    
    C32(1,:) = 0; C32(:,1) = 0;
    C32(N1/2+1,:) = 0; C32(:,N2/2+1) = 0;
    
    C33(1,:) = 0; C33(:,1) = 0;
    C33(N1/2+1,:) = 0; C33(:,N2/2+1) = 0;
    
    h5write(h5C, '/C11', C11, [1 1 k], [N1 N2 1]);
    h5write(h5C, '/C21', C21, [1 1 k], [N1 N2 1]);
    h5write(h5C, '/C31', C31, [1 1 k], [N1 N2 1]);
    h5write(h5C, '/C22', C22, [1 1 k], [N1 N2 1]);
    h5write(h5C, '/C32', C32, [1 1 k], [N1 N2 1]);
    h5write(h5C, '/C33', C33, [1 1 k], [N1 N2 1]);
    
    fprintf('C finished for k=%u.\n', k);
end

%% FFT
h5create(h5temp, '/u_re', [N1 N2 N3]);
h5create(h5temp, '/u_im', [N1 N2 N3]);
h5create(h5temp, '/v_re', [N1 N2 N3]);
h5create(h5temp, '/v_im', [N1 N2 N3]);
h5create(h5temp, '/w_re', [N1 N2 N3]);
h5create(h5temp, '/w_im', [N1 N2 N3]);

for k = 1:N3
    C11 = h5read(h5C, '/C11', [1 1 k], [N1 N2 1]);
    C21 = h5read(h5C, '/C21', [1 1 k], [N1 N2 1]);
    C31 = h5read(h5C, '/C31', [1 1 k], [N1 N2 1]);
    C22 = h5read(h5C, '/C22', [1 1 k], [N1 N2 1]);
    C32 = h5read(h5C, '/C32', [1 1 k], [N1 N2 1]);
    C33 = h5read(h5C, '/C33', [1 1 k], [N1 N2 1]);

    Phiu = h5read(h5temp, '/Phiu', [1 1 k], [N1 N2 1]);
    Phiv = h5read(h5temp, '/Phiv', [1 1 k], [N1 N2 1]);
    Phiw = h5read(h5temp, '/Phiw', [1 1 k], [N1 N2 1]);

    Au = C11 .* exp(1i * Phiu);
    Av = C21 .* exp(1i * Phiu) + C22 .* exp(1i * Phiv);
    Aw = C31 .* exp(1i * Phiu) + C32 .* exp(1i * Phiv) + C33 .* exp(1i * Phiw);

    u = N1 * N2 * ifft(ifft(Au, [], 1), [], 2);
    v = N1 * N2 * ifft(ifft(Av, [], 1), [], 2);
    w = N1 * N2 * ifft(ifft(Aw, [], 1), [], 2);

    h5write(h5temp, '/u_re', real(u), [1 1 k], [N1 N2 1]);
    h5write(h5temp, '/u_im', imag(u), [1 1 k], [N1 N2 1]);
    h5write(h5temp, '/v_re', real(v), [1 1 k], [N1 N2 1]);
    h5write(h5temp, '/v_im', imag(v), [1 1 k], [N1 N2 1]);
    h5write(h5temp, '/w_re', real(w), [1 1 k], [N1 N2 1]);
    h5write(h5temp, '/w_im', imag(w), [1 1 k], [N1 N2 1]);
    
    fprintf('FFT along k1 and k2 finished for k3 index = %u\n', k);
end

h5file = 'windField.h5';

if ~exist(h5file, 'file')
    h5create(h5file, '/u', [N1 N2 N3]);
    h5create(h5file, '/v', [N1 N2 N3]);
    h5create(h5file, '/w', [N1 N3 N3]);
end

for j = 1:N2
    u_re = h5read(h5temp, '/u_re', [1 j 1], [N1 1 N3]);
    u_im = h5read(h5temp, '/u_im', [1 j 1], [N1 1 N3]);
    v_re = h5read(h5temp, '/v_re', [1 j 1], [N1 1 N3]);
    v_im = h5read(h5temp, '/v_im', [1 j 1], [N1 1 N3]);
    w_re = h5read(h5temp, '/w_re', [1 j 1], [N1 1 N3]);
    w_im = h5read(h5temp, '/w_im', [1 j 1], [N1 1 N3]);
    
    u = N3 * ifft(u_re + 1i * u_im, [], 3);
    v = N3 * ifft(v_re + 1i * v_im, [], 3);
    w = N3 * ifft(w_re + 1i * w_im, [], 3);
    
    fprintf('FFT along k3 finished for k2 index = %u\n', j);
    
    h5write(h5file, '/u', real(u), [1 j 1], [N1 1 N3]);
    h5write(h5file, '/v', real(v), [1 j 1], [N1 1 N3]);
    h5write(h5file, '/w', real(w), [1 j 1], [N1 1 N3]);
end

h5writeatt(h5file, '/', 'gamma', gamma);
h5writeatt(h5file, '/', 'L', L);
h5writeatt(h5file, '/', 'ae23', ae23);

h5writeatt(h5file, '/', 'Lx10', Lx10);
h5writeatt(h5file, '/', 'Lx20', Lx20);
h5writeatt(h5file, '/', 'Lx30', Lx30);

h5writeatt(h5file, '/', 'N1', N1);
h5writeatt(h5file, '/', 'N2', N2);
h5writeatt(h5file, '/', 'N3', N3);
