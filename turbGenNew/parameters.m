%% Coordinates system for CFD computational domain
% x-axis: main flow direction (along-wind, longitudinal)
% y-axis: transverse (cross-wind)
% z-axis: vertical (z coordinates represent height above ground)

%% Mean velocity
% h0u:    Reference height for the mean velocity
% Uh:     Mean velocity at h0u
% alphau: Power low exponent of the mean velocity
h0u = 0.364;
% alphau = 0.326;
alphau = 0;
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
nt = 20481;
Td = nt * dt;

%% Frequency segments
% N:      Number of frequency segments
% df:     frequency step
df = 1 / Td;
N = (nt+1)/2;
fmax = (N-1)*df;

%% Wavenumber segments
% M:      Number of wavenumber segments
% kmax:   Maximum wavenumber
M = 3000;
kmax = 3000;

%% Sample test grid
% dh = 0.001;
% dr = 0.2;
% x0 = 0;
% y0 = 0;
% z0 = 0.1;
% GRID = [x0 y0 z0;        % Pt. 1
%         %x0 y0 z0+2*dh;  % Pt. 2
%         %x0-dh y0 z0+dh; % Pt. 3
%         %x0+dh y0 z0+dh; % Pt. 4
%         %x0 y0-dh z0+dh; % Pt. 5
%         %x0 y0+dh z0+dh; % Pt. 6
%         %x0 y0 z0+dh;    % Pt. 7
%         x0 y0 z0+dr;     % Pt. 8
%     ];
GRID = [zeros(10,1) zeros(10,1) (0.05:0.1:1)'];