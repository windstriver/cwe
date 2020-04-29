clear; clc;
%% Spectral tensor parameter
% Kaimal model
gamma = 3.9;
z = 40; %[m]
L = 0.59 * z; %[m]
ae23 = 0.86; %[m^(4/3) s^(-1)

u_tau = sqrt(ae23 * z.^(2/3) / 3.2);

%% Wavenumber grid and spectral tensor grid
k1u = 8*pi/L;
k2u = 8*pi/L;
k3u = 8*pi/L;

N1 = 513;
N2 = 513;
N3 = 513;

k1_vec = linspace(-k1u, k1u, N1)';
k2_vec = linspace(-k2u, k2u, N2)';
k3_vec = linspace(-k3u, k3u, N3)';

delta_k1 = 2*k1u / (N1-1);
delta_k2 = 2*k2u / (N2-1);
delta_k3 = 2*k3u / (N3-1);

h5file = 'Phi.h5';
if ~exist(h5file, 'file')
    h5create(h5file, '/Phi11', [N1 N2 N3]);
    h5create(h5file, '/Phi22', [N1 N2 N3]);
    h5create(h5file, '/Phi33', [N1 N2 N3]);
    h5create(h5file, '/Phi12', [N1 N2 N3]);
    h5create(h5file, '/Phi13', [N1 N2 N3]);
    h5create(h5file, '/Phi23', [N1 N2 N3]);
end

for k = 1:N3
    [k1_grid, k2_grid, k3_grid] = ndgrid(k1_vec, k2_vec, k3_vec(k));

    [Phi11, Phi22, Phi33, Phi12, Phi13, Phi23] = ...
        spectralTensor(k1_grid, k2_grid, k3_grid, gamma, L, ae23);
    
    % Save spectral tensor to HDF5 database
    h5write(h5file, '/Phi11', Phi11, [1 1 k], [N1 N2 1]);
    h5write(h5file, '/Phi22', Phi22, [1 1 k], [N1 N2 1]);
    h5write(h5file, '/Phi33', Phi33, [1 1 k], [N1 N2 1]);
    h5write(h5file, '/Phi12', Phi12, [1 1 k], [N1 N2 1]);
    h5write(h5file, '/Phi13', Phi13, [1 1 k], [N1 N2 1]);
    h5write(h5file, '/Phi23', Phi23, [1 1 k], [N1 N2 1]);
    
    fprintf('Spectral tensor for k3 index = %u finished.\n', k);
end

h5writeatt(h5file, '/', 'gamma', gamma);
h5writeatt(h5file, '/', 'L', L);
h5writeatt(h5file, '/', 'ae23', ae23);

h5writeatt(h5file, '/', 'N1', N1);
h5writeatt(h5file, '/', 'N2', N2);
h5writeatt(h5file, '/', 'N3', N3);

h5writeatt(h5file, '/', 'delta_k1', delta_k1);
h5writeatt(h5file, '/', 'delta_k2', delta_k2);
h5writeatt(h5file, '/', 'delta_k3', delta_k3);
h5writeatt(h5file, '/', 'delta_k1', delta_k1);

%% Comparison of 1D spectrum
% F1 = trapz(trapz(Phi11, 2) * delta_k2, 3) * delta_k3;
% F2 = trapz(trapz(Phi22, 2) * delta_k2, 3) * delta_k3;
% F3 = trapz(trapz(Phi33, 2) * delta_k2, 3) * delta_k3;
% 
% k1_vec_pos = k1_vec(k1_vec > 0);
% F1_pos = F1(k1_vec > 0);
% F2_pos = F2(k1_vec > 0);
% F3_pos = F3(k1_vec > 0);
% 
% k1_vec_pos_norm = k1_vec_pos * z / (2*pi);
% 
% % k1_vec_pos_norm = 1e-3:1e-3:10;
% 
% figure;
% % longitudinal
% loglog(k1_vec_pos_norm*2*pi, 52.5*k1_vec_pos_norm .* ((1+33*k1_vec_pos_norm).^(-5/3)),...
%     'r-', 'DisplayName', 'Kaimal F1');
% hold on;
% loglog(k1_vec_pos_norm*2*pi, k1_vec_pos.*F1_pos/u_tau^2, ...
%     'r-', 'Marker', '.', 'DisplayName', 'Mann F1');
% % lateral
% loglog(k1_vec_pos_norm*2*pi, 8.5*k1_vec_pos_norm .* ((1+9.5*k1_vec_pos_norm).^(-5/3)),...
%     'g--', 'DisplayName', 'Kaimal F2');
% loglog(k1_vec_pos_norm*2*pi, k1_vec_pos.*F2_pos/u_tau^2, ...
%     'g--', 'Marker', 's', 'DisplayName', 'Mann F2');
% % lateral
% loglog(k1_vec_pos_norm*2*pi, 1.05*k1_vec_pos_norm ./ (1+5.3*k1_vec_pos_norm.^(5/3)),...
%     'b-.', 'DisplayName', 'Kaimal F3');
% loglog(k1_vec_pos_norm*2*pi, k1_vec_pos.*F3_pos/u_tau^2, ...
%     'b-.', 'Marker', '^', 'DisplayName', 'Mann F3');
% hold off;
% grid on;
% legend();
% xlabel('$ k_1 z $', 'Interpreter', 'latex');
% ylabel('$ k_1 F(k_1) / u_{\tau}^2 $', 'Interpreter', 'latex');
% xticks([0.1, 1, 10, 100]);
% yticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5]);

% Cross spectra of u and w
% cov_12 = trapz(trapz(trapz(Phi12, 1), 2), 3) * delta_k1 * delta_k2 * delta_k3;
% cov_23 = trapz(trapz(trapz(Phi23, 1), 2), 3) * delta_k1 * delta_k2 * delta_k3;
% cov_13 = trapz(trapz(trapz(Phi13, 1), 2), 3) * delta_k1 * delta_k2 * delta_k3;
% 
% F13 = trapz(trapz(Phi13, 2), 3) * delta_k2 * delta_k3;
% F13_pos = F13(k1_vec > 0);
% 
% semilogx(k1_vec_pos_norm*2*pi, (-7)*k1_vec_pos_norm .* ((1+9.6*k1_vec_pos_norm).^(-7/3)),...
%     'k-', 'DisplayName', 'Kaimal F13');
% hold on;
% semilogx(k1_vec_pos_norm*2*pi, k1_vec_pos.*F13_pos/u_tau^2, ...
%     'r-', 'Marker', '.', 'DisplayName', 'Mann F13');
% hold off;
% grid on;
% legend();
% xlabel('$ k_1 z $', 'Interpreter', 'latex');
% ylabel('$ k_1 F_{13}(k_1) / u_{\tau}^2 $', 'Interpreter', 'latex');

% Coherency
% Coh13_pos = abs(F13_pos) ./ sqrt(F1_pos) ./ sqrt(F3_pos);
% semilogx(k1_vec_pos, Coh13_pos)
