%% Coordinates system for CFD computational domain
% x-axis: streamwise
% y-axis: wall-normal
% z-axis: spanwise
% Units: m-sec

%% Read grid of the inflow plane from OpenFOAM
GRID = csvread('../testCase/constant/meshInfo/faceCentresInlet.csv');
nd = size(GRID,1);  % overall number of points

%% Extract the coordinates
% Y coordinate is the height above ground
X = GRID(:,2);
Y = GRID(:,3);
Z = GRID(:,4);

%% Mean velocity profile
% H_ref:    Reference height for the mean velocity
% U_ref:    Mean velocity at H_ref
% y0:       Roughness length in model scale
% alpha:    Exposure exponent
H_ref = 0.5;  %[m]
U_ref = 11.1438;     %[m/s]
% y0 = 0.02/300;  %[m]
alpha = 1/4;
% Uav = U_ref / log(H_ref/y0) * log(Y/y0);
Uav = U_ref * (Y/H_ref).^alpha;

%% Turbulent intensity profile
% Iu = 1.265 * 1.00 ./ log(Y/y0);
% Iv = 1.075 * 0.08 *  ones(nd,1);
% Iw = 1.060 * 0.88 ./ log(Y/(2.5e-5));
Iu = 1.10 * 0.116 * (Y/H_ref).^(-alpha-0.05);
Iv = 1.00 * 0.550 * Iu;
Iw = 1.00 * 0.780 * Iu;

%% Integrale length scale profile
% Lu = 2.0 * 2.1982*H_ref*(Y/H_ref).^(0.1209);
% Lv = 0.7137*H_ref*(Y/H_ref).^(0.3652);
% Lw = 0.6229*H_ref*(Y/H_ref).^(0.3505);
Lu = 1.65 * 1.480 * H_ref * (Y/H_ref).^(0.4);
Lv = 0.58 * 0.082 * Lu;
Lw = 0.67 * 0.237 * Lu;

%% Coherency decay constants
% Cxyz:   Coherency decay constants in x, y and z directions [1 3] matrix
Cxyz = [10 10 10];

%% Simulate time steps
% dt:    time step
% nt:    number of time steps
% Td:    total simulated time
dt = 5e-4;
nt = 120001;
Td = nt * dt;

%% Frequency segments
% N:      number of frequency segments
% df:     frequency step
df = 1 / Td;
N = (nt+1) / 2;
fmax = (N-1) * df;

%% Wavenumber segments
% M:      Number of wavenumber segments
% kmax:   Maximum wavenumber
M = 1000;
kmax = 1000;

%% Extract the frequency and time vector
fvec = (0:(N-1)) * df;   % frequency vector
tvec = dt * (0:(nt-1))'; % time vector

%% Extract the wavenumber vector
kmin = kmax / M;          % min wavenumber
dk = (kmax-kmin) / (M-1); % wavenumber step
kvec = (kmin:dk:kmax)';   % wavenumber vector

%% Calculate the von Karmon spectrum matrices
% row index is for frequency segments, {0, df, ..., (N-1)*df}
% column index is for different points, [P1, P2, ..., Pnd]
Su0 = 4*ones(N,1)*((Iu'.*Uav').^2.*(Lu'./Uav')) ./ ...
    (1 + 70.8*(fvec'*(Lu'./Uav')).^2).^(5/6);
Sv0 = 4*ones(N,1)*((Iv'.*Uav').^2.*(Lv'./Uav')) .* ...
    (1 + 188.4*(2*fvec'*(Lv'./Uav')).^2) ./ ...
    (1 + 70.8*(2*fvec'*(Lv'./Uav')).^2).^(11/6);
Sw0 = 4*ones(N,1)*((Iw'.*Uav').^2.*(Lw'./Uav')) .* ...
    (1 + 188.4*(2*fvec'*(Lw'./Uav')).^2) ./ ...
    (1 + 70.8*(2*fvec'*(Lw'./Uav')).^2).^(11/6);

% Create the HDF5 file
hdf5File = 'inflowTurb2.h5';
if exist(hdf5File, 'file') == 0
    run('hdf5Creation.m');
end

