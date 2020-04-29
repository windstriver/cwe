function [Phi11, Phi22, Phi33, Phi12, Phi13, Phi23] = ...
    spectralTensor(k1_grid, k2_grid, k3_grid, gamma, L, ae23)

%SPECTRALTENSOR Spectral tensor model of Mann.
% Input:
% ---- k1_grid, wavenumber component 1 after distorsion by shear
% ---- k2_grid, wavenumber component 2 after distorsion by shear
% ---- k3_grid, wavenumber component 3 after distorsion by shear
% ---- gamma,   parameter of the sheared spectral tensor
% ---- L    ,   length scale of the spectral velocity tensor [m]
% ---- ae23 ,   \alpha \epsilon^(2/3)
%               \alpha,   3D Kolmogorov constant (about 1.7)
%               \epsilon, dissipation of turbulent kinetic energy
% Output:
% ---- Phi11, spectral velocity tensor component (1,1)
% ---- Phi22, spectral velocity tensor component (2,2)
% ---- Phi33, spectral velocity tensor component (3,3)
% ---- Phi12, spectral velocity tensor component (1,2)
% ---- Phi13, spectral velocity tensor component (1,3)
% ---- Phi23, spectral velocity tensor component (2,3)
%             1,2,3 denotes longitudinal, lateral, and vertical directions

VSMALL = 1e-6;

mask_k1_0 = (k1_grid == 0);
k1_grid(mask_k1_0) = VSMALL;

mask_k2_0 = (k2_grid == 0);
k2_grid(mask_k2_0) = VSMALL;

mask_k3_0 = (k3_grid == 0);
k3_grid(mask_k3_0) = VSMALL;

kmag_grid = sqrt(k1_grid.^2+k2_grid.^2+k3_grid.^2);

beta_grid = gamma ./ (kmag_grid .* L).^(2/3) ./ ...
    sqrt(hypergeom([1/3 17/6], 4/3, -(kmag_grid .* L).^(-2)));

% wavenumber component 3 before distortion
k30_grid = k3_grid + beta_grid .* k1_grid;
k0mag_grid = sqrt(k1_grid.^2+k2_grid.^2+k30_grid.^2);

C1_grid = beta_grid .* k1_grid .^2 .* ...
    (k0mag_grid .^2 - 2 * k30_grid .^2 + beta_grid .* k1_grid .* k30_grid) ./ ...
    (kmag_grid .^2 .* (k1_grid .^2 + k2_grid .^2));

C2_grid = k2_grid .* k0mag_grid .^2 ./ ...
    (k1_grid .^2 + k2_grid .^2) .^(3/2) .* ...
    atan(beta_grid .* k1_grid .* sqrt(k1_grid .^2 + k2_grid .^2) ./ ...
    (k0mag_grid .^2 - k30_grid .* k1_grid .* beta_grid));

zeta1_grid = C1_grid - k2_grid ./ k1_grid .* C2_grid;
zeta2_grid = k2_grid ./ k1_grid .* C1_grid + C2_grid;

E_k0_grid = ae23 * L^(5/3) .* (L*k0mag_grid).^4 ./ ...
    (1+(L*k0mag_grid) .^2) .^(17/6);

Phi11 = E_k0_grid ./ (4*pi*k0mag_grid .^4) .* ...
    (k0mag_grid.^2 - k1_grid.^2 - 2*k1_grid.*k30_grid.*zeta1_grid + ...
    (k1_grid.^2 + k2_grid.^2) .* zeta1_grid.^2);

Phi22 = E_k0_grid ./ (4*pi*k0mag_grid .^4) .* ...
    (k0mag_grid.^2 - k2_grid.^2 - 2*k2_grid.*k30_grid.*zeta2_grid + ...
    (k1_grid.^2 + k2_grid.^2) .* zeta2_grid.^2);

Phi33 = E_k0_grid ./ (4*pi*kmag_grid .^4) .* ...
    (k1_grid.^2 + k2_grid.^2);

Phi12 = E_k0_grid ./ (4*pi*k0mag_grid .^4) .* ...
    (-k1_grid .* k2_grid - k1_grid .* k30_grid .* zeta2_grid - ...
    k2_grid .* k30_grid .* zeta1_grid + ...
    (k1_grid .^2 + k2_grid .^2) .* zeta1_grid .* zeta2_grid);

Phi13 = E_k0_grid ./ (4*pi*k0mag_grid .^2 .* kmag_grid .^2) .* ...
    (-k1_grid .* k30_grid + ...
    (k1_grid .^ 2 + k2_grid .^ 2) .* zeta1_grid);

Phi23 = E_k0_grid ./ (4*pi*k0mag_grid .^2 .* kmag_grid .^2) .* ...
    (-k2_grid .* k30_grid + ...
    (k1_grid .^ 2 + k2_grid .^ 2) .* zeta2_grid);

Phi11(mask_k1_0) = 0;
Phi11(mask_k2_0) = 0;
Phi11(mask_k3_0) = 0;

Phi22(mask_k1_0) = 0;
Phi22(mask_k2_0) = 0;
Phi22(mask_k3_0) = 0;

Phi33(mask_k1_0) = 0;
Phi33(mask_k2_0) = 0;
Phi33(mask_k3_0) = 0;

Phi12(mask_k1_0) = 0;
Phi12(mask_k2_0) = 0;
Phi12(mask_k3_0) = 0;

Phi13(mask_k1_0) = 0;
Phi13(mask_k2_0) = 0;
Phi13(mask_k3_0) = 0;

Phi23(mask_k1_0) = 0;
Phi23(mask_k2_0) = 0;
Phi23(mask_k3_0) = 0;

end

