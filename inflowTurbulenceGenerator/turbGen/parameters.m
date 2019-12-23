%% Coordinates system
% x-axis: streamwise
% y-axis: wall-normal
% z-axis: spanwise
% Units: m-sec

%% Parameters
% H_ref:    Reference height for the mean velocity
% U_ref:    Mean velocity at H_ref
H_ref = 0.5; %[m]
U_ref = 11.1438; %[m/s]

%% Grid used to generate the inflow
% GRID = [ (6:10)', -5*H_ref*ones(5,1), (6:10)'*0.1*H_ref, zeros(5,1)];
GRID = csvread('../testCase/constant/meshInfo/faceCentresinlet.csv');
nd = size(GRID,1);  % overall number of points

%% Extract the coordinates
% Y coordinate is the vertical height from the wall
X = GRID(:,1); Y = GRID(:,2); Z = GRID(:,3);

%% Mean velocity profile from TPU database
alpha = 1/4;
Uav = U_ref * (Y/H_ref).^alpha;

%% Turbulent intensity profile
Iu = 0.116 * (Y/H_ref).^(-alpha-0.05);
Iv = 0.55 * Iu;
Iw = 0.78 * Iu;

%% Integrale length scale profile for TPU database
Lu = 1.0 * H_ref * ones(nd,1);
Lv = 0.082 * Lu;
Lw = 0.237 * Lu;

%% Coherency decay constants
% Cxyz:   Coherency decay constants in x, y and z directions [1 3] matrix
Cxyz = [10 10 10];

%% Simulate time steps
% dt:    time step
% nt:    number of time steps
% Td:    total simulated time
dt = 5e-4;
nt = 48001;
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
M = 5000;
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
hdf5File = 'inflowTurb_samp_freq_2e3.h5';
try
    run('hdf5Creation.m');
catch
end

