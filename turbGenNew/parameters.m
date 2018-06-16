%% Coordinates system for CFD computational domain
% x-axis: main flow direction (along-wind, longitudinal)
% y-axis: transverse (cross-wind)
% z-axis: vertical (z coordinates represent height above ground)
% Units: m-sec

%% Sample test grid
GRID = [(0:9)' zeros(10,1) zeros(10,1) (0.05:0.1:1)'];
nd = size(GRID,1);  % overall number of points

%% Mean velocity
% h0u:    Reference height for the mean velocity
% Uh:     Mean velocity at h0u
% alphau: Power low exponent of the mean velocity
h0u = 0.364;
alphau = 0.326;
% alphau = 0;
Uh = 10.0;

%% Turbulent intensity
% h0I:    Reference height for the turbulent intensity
% Iuh:    Longitudinal turbulent intensity at h0I
% Ivh:    Transverse turbulent intensity at h0I
% Iwh:    Vertical turbulent intensity at h0I
% dIu:    Power low exponent of the longitudinal turbulent intensity
% dIv:    Power low exponent of the transverse turbulent intensity
% dIw:    Power low exponent of the vertical turbulent intensity
h0I = 0.364;
Iuh = 0.208;
Ivh = 0.182;
Iwh = 0.152;
% dIu = -0.191;
dIu = 0;
% dIv = -0.123;
dIv = 0;
% dIw = -0.005;
dIw = 0;

%% Turbulence length scale
% h0L:    Reference height for the length scale
% Luh:    Longitudinal length scale at h0L
% Lvh:    Transverse length scale at h0L
% Lwh:    Vertical length scale at h0L
% dLu:    Power low exponent of the longitudinal length scale
% dLv:    Power low exponent of the transverse length scale
% dLw:    Power low exponent of the vertical length scale
h0L = 0.254;
Luh = 0.302;
Lvh = 0.0815;
Lwh = 0.0326;
% dLu = 0.473;
dLu = 0;
% dLv = 0.881;
dLv = 0;
% dLw = 1.539;
dLw = 0;

%% Coherency decay constants
% Cxyz:   Coherency decay constants in x, y and z directions [1 3] matrix
Cxyz = [10 10 10];

%% Simulate time steps
% dt:    time step
% nt:    number of time steps
% Td:    total simulated time
dt = 0.001;
nt = 40961;
Td = nt * dt;

%% Frequency segments
% N:      Number of frequency segments
% df:     frequency step
df = 1 / Td;
N = (nt+1) / 2;
fmax = (N-1) * df;

%% Wavenumber segments
% M:      Number of wavenumber segments
% kmax:   Maximum wavenumber
M = 3000;
kmax = 3000;

%% Extract the coordinates
% inflow plane is z-plane
X = GRID(:,2);
Y = GRID(:,3);
Z = GRID(:,4);

%% Extract the frequency and time vector
fvec = (0:(N-1)) * df; % frequency vector
tvec = dt * (0:(nt-1))'; % time vector

%% Extract the wavenumber vector
kmin = kmax / M; % min wavenumber
dk = (kmax-kmin) / (M-1); % wavenumber step
kvec = (kmin:dk:kmax)'; % wavenumber vector

%% Calculate average velocity, turbulent Intensity, length scale profiles
% Y coordinate is the height above ground
Uav = Uh*(Z/h0u).^alphau;
Iu = Iuh*(Z/h0I).^dIu;
Iv = Ivh*(Z/h0I).^dIv;
Iw = Iwh*(Z/h0I).^dIw;
Lu = Luh*(Z/h0L).^dLu;
Lv = Lvh*(Z/h0L).^dLv;
Lw = Lwh*(Z/h0L).^dLw;

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

%% Create HDF5 database to store simulated data
% HDF5 library uses C-style ordering for multidimensional arrays
% MATLAB uses FORTRAN-style ordering
% Dataset size in HDF5 database for use in OpenFOAM
% GRID: [nd 4]    PtNum, x, y, z
% TIME: [nt]
% UMEAN: [nd]
% U,V,W: [nd nt]
hdf5File = 'lesinlet.h5';
if (exist(hdf5File,'file') == 2)
    delete(hdf5File);
end

h5create(hdf5File,'/GRID',fliplr([nd 4]));
h5create(hdf5File,'/TIME',nt);
h5create(hdf5File,'/UMEAN',nd);
h5create(hdf5File,'/U',fliplr([nd nt]));
h5create(hdf5File,'/V',fliplr([nd nt]));
h5create(hdf5File,'/W',fliplr([nd nt]));
