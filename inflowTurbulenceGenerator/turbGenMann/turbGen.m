clear; clc;
%% Spectral tensor parameter
% Kaimal model
gamma = 3.9;
z = 40; %[m]
L = 0.59 * z; %[m]
ae23 = 0.86; %[m^(4/3) s^(-1)

u_tau = sqrt(ae23 * z.^(2/3) / 3.2);

%% Simulation box parameters
Lx10 =  32 * L;
Lx20 =  16 * L;
Lx30 =  16 * L;

N1 = 1024;
N2 =  512;
N3 =  512;

delta_x1 = Lx10 / N1;
delta_x2 = Lx20 / N2;
delta_x3 = Lx30 / N3;

delta_k1 = 2*pi / Lx10;
delta_k2 = 2*pi / Lx20;
delta_k3 = 2*pi / Lx30;

%% Random phases
Phiu = zeros(N1, N2, N3);
Phiv = zeros(N1, N2, N3);
Phiw = zeros(N1, N2, N3);

Phiu(2:N1/2,2:N2/2,2:N3/2) = 2*pi*rand(N1/2-1, N2/2-1, N3/2-1);
Phiu(N1/2+2:N1,N2/2+2:N2,N3/2+2:N3) = -Phiu(N1/2:-1:2,N2/2:-1:2,N3/2:-1:2);

Phiu(2:N1/2,N2/2+2:N2,2:N3/2) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiu(N1/2+2:N1,2:N2/2,N3/2+2:N3) = -Phiu(N1/2:-1:2,N2:-1:N2/2+2,N3/2:-1:2);

Phiu(2:N1/2,2:N2/2,N3/2+2:N3) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiu(N1/2+2:N1,N2/2+2:N2,2:N3/2) = -Phiu(N1/2:-1:2,N2/2:-1:2,N3:-1:N3/2+2);

Phiu(2:N1/2,N2/2+2:N2,N3/2+2:N3) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiu(N1/2+2:N1,2:N2/2,2:N3/2) = -Phiu(N1/2:-1:2,N2:-1:N2/2+2,N3:-1:N3/2+2);


Phiv(2:N1/2,2:N2/2,2:N3/2) = 2*pi*rand(N1/2-1, N2/2-1, N3/2-1);
Phiv(N1/2+2:N1,N2/2+2:N2,N3/2+2:N3) = -Phiv(N1/2:-1:2,N2/2:-1:2,N3/2:-1:2);

Phiv(2:N1/2,N2/2+2:N2,2:N3/2) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiv(N1/2+2:N1,2:N2/2,N3/2+2:N3) = -Phiv(N1/2:-1:2,N2:-1:N2/2+2,N3/2:-1:2);

Phiv(2:N1/2,2:N2/2,N3/2+2:N3) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiv(N1/2+2:N1,N2/2+2:N2,2:N3/2) = -Phiv(N1/2:-1:2,N2/2:-1:2,N3:-1:N3/2+2);

Phiv(2:N1/2,N2/2+2:N2,N3/2+2:N3) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiv(N1/2+2:N1,2:N2/2,2:N3/2) = -Phiv(N1/2:-1:2,N2:-1:N2/2+2,N3:-1:N3/2+2);


Phiw(2:N1/2,2:N2/2,2:N3/2) = 2*pi*rand(N1/2-1, N2/2-1, N3/2-1);
Phiw(N1/2+2:N1,N2/2+2:N2,N3/2+2:N3) = -Phiw(N1/2:-1:2,N2/2:-1:2,N3/2:-1:2);

Phiw(2:N1/2,N2/2+2:N2,2:N3/2) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiw(N1/2+2:N1,2:N2/2,N3/2+2:N3) = -Phiw(N1/2:-1:2,N2:-1:N2/2+2,N3/2:-1:2);

Phiw(2:N1/2,2:N2/2,N3/2+2:N3) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiw(N1/2+2:N1,N2/2+2:N2,2:N3/2) = -Phiw(N1/2:-1:2,N2/2:-1:2,N3:-1:N3/2+2);

Phiw(2:N1/2,N2/2+2:N2,N3/2+2:N3) = 2*pi*rand(N1/2-1,N2/2-1,N3/2-1);
Phiw(N1/2+2:N1,2:N2/2,2:N3/2) = -Phiw(N1/2:-1:2,N2:-1:N2/2+2,N3:-1:N3/2+2);

%% Preparation of Fourier transform of velocity vector random process
k1_vec = [0, 1:(N1/2), -(N1/2-1):(-1)] * delta_k1;
k2_vec = [0, 1:(N2/2), -(N2/2-1):(-1)] * delta_k2;
k3_vec = [0, 1:(N3/2), -(N3/2-1):(-1)] * delta_k3;

[k1_grid, k2_grid, k3_grid] = ndgrid(k1_vec, k2_vec, k3_vec);

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
C33 = sqrt(RHS33 - C32 .^2 - C31 .^2);

C11(1,:,:) = 0; C11(:,1,:) = 0; C11(:,:,1) = 0;
C11(N1/2+1,:,:) = 0; C11(:,N2/2+1,:) = 0; C11(:,:,N3/2+1) = 0;

C22(1,:,:) = 0; C22(:,1,:) = 0; C22(:,:,1) = 0;
C22(N1/2+1,:,:) = 0; C22(:,N2/2+1,:) = 0; C22(:,:,N3/2+1) = 0;

C33(1,:,:) = 0; C33(:,1,:) = 0; C33(:,:,1) = 0;
C33(N1/2+1,:,:) = 0; C33(:,N2/2+1,:) = 0; C33(:,:,N3/2+1) = 0;

C21(1,:,:) = 0; C21(:,1,:) = 0; C21(:,:,1) = 0;
C21(N1/2+1,:,:) = 0; C21(:,N2/2+1,:) = 0; C21(:,:,N3/2+1) = 0;

C31(1,:,:) = 0; C31(:,1,:) = 0; C31(:,:,1) = 0;
C31(N1/2+1,:,:) = 0; C31(:,N2/2+1,:) = 0; C31(:,:,N3/2+1) = 0;

C32(1,:,:) = 0; C32(:,1,:) = 0; C32(:,:,1) = 0;
C32(N1/2+1,:,:) = 0; C32(:,N2/2+1,:) = 0; C32(:,:,N3/2+1) = 0;

% independent Gaussian random variables
Au = C11 .* exp(1i * Phiu);
Av = C21 .* exp(1i * Phiu) + C22 .* exp(1i * Phiv);
Aw = C31 .* exp(1i * Phiu) + C32 .* exp(1i * Phiv) + C33 .* exp(1i * Phiw);

%% FFT
u = real(N1*N2*N3*ifftn(Au));
v = real(N1*N2*N3*ifftn(Av));
w = real(N1*N2*N3*ifftn(Aw));

h5file = 'windField.h5';
if ~exist(h5file, 'file')
    h5create(h5file, '/u', size(u));
    h5create(h5file, '/v', size(v));
    h5create(h5file, '/w', size(w));
end

h5write(h5file, '/u', u);
h5write(h5file, '/v', v);
h5write(h5file, '/w', w);

h5writeatt(h5file, '/', 'gamma', gamma);
h5writeatt(h5file, '/', 'L', L);
h5writeatt(h5file, '/', 'ae23', ae23);

h5writeatt(h5file, '/', 'Lx10', Lx10);
h5writeatt(h5file, '/', 'Lx20', Lx20);
h5writeatt(h5file, '/', 'Lx30', Lx30);

h5writeatt(h5file, '/', 'N1', N1);
h5writeatt(h5file, '/', 'N2', N2);
h5writeatt(h5file, '/', 'N3', N3);
