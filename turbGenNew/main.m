%% Run parameters.m
parameters

%% Extract the coordinates
% inflow plane is z-plane
X = GRID(:,2);
Y = GRID(:,3);
Z = GRID(:,4);

%% Extract the frequency and time vector
fvec = (0:(N-1))*df; % frequency vector
tvec = dt * (0:(nt-1))'; % time vector

%% Extract the wavenumber vector
kmin = kmax/M; % min wavenumber
dk = (kmax-kmin)/(M-1); % wavenumber step
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

%% Generate random phase matrices
PSI = 2*pi*rand(M,N);
PHI = 2*pi*rand(M,N);

%% Start generating turbulence
% Initialization
% Velocity time series
% row: time step index
% column: points index
u = zeros(nt,nd);
v = zeros(nt,nd);
w = zeros(nt,nd);

parfor i = 1:nd    % i: points index
    % Wavenumber-Freq Spectrum matrix at a single point
    Sukf = zeros(M,N);    % Suu(km, fn)
    Svkf = zeros(M,N);    % Svv(km, fn)
    Swkf = zeros(M,N);    % Sww(km, fn)

    % Amplitude for wave(km, fn)
    Umn = zeros(M,N);
    Vmn = zeros(M,N);
    Wmn = zeros(M,N);

    % kxmn for divergence-free condition
    kxmn = zeros(M,N);

    % Phase for wave(fn)
    Phmn = zeros(M,N);

    % FFT matrix for velocity time series u, v, w
    % row: frequency index
    % 0, df, 2*df, ... , (N-1)*df
    % N*df, (N+1)*df, ... , (2*N-2)*df
    % column: wavenumber index
    Uf = zeros(2*N-1, M);
    Vf = zeros(2*N-1, M);
    Wf = zeros(2*N-1, M);

    % Wave-freq. Spectrum matrix for point i
    Sukf = 2/pi * ones(M,1) * (Su0(:,i)'.*(Cxyz(1)*fvec/Uav(i))) ./ ...
        ((Cxyz(1)*fvec/Uav(i)).^2 + kvec.^2);
    Svkf = 2/pi * ones(M,1) * (Sv0(:,i)'.*(Cxyz(2)*fvec/Uav(i))) ./ ...
        ((Cxyz(2)*fvec/Uav(i)).^2 + kvec.^2);
    Swkf = 2/pi * ones(M,1) * (Sw0(:,i)'.*(Cxyz(3)*fvec/Uav(i))) ./ ...
        ((Cxyz(3)*fvec/Uav(i)).^2 + kvec.^2);
    % Compare std. calculated by integration of wave-freq. spectrum with
    % theory
    fprintf('Pt. %3d (%6.3f, %6.3f, %6.3f):\n', i, X(i), Y(i), Z(i));
    fprintf('Std. of along-wind speed by theory: %6.3f\n', Iu(i)*Uav(i));
    fprintf('Std. of along-wind speed by wave-freq spectrum: %6.3f\n',...
        sqrt(sum(sum(Sukf))*df*dk));
    fprintf('Std. of cross-wind speed by theory: %6.3f\n', Iv(i)*Uav(i));
    fprintf('Std. of cross-wind speed by wave-freq spectrum: %6.3f\n',...
        sqrt(sum(sum(Svkf))*df*dk));
    fprintf('Std. of vertical wind speed by theory: %6.3f\n', Iw(i)*Uav(i));
    fprintf('Std. of vertical wind speed by wave-freq spectrum: %6.3f\n',...
        sqrt(sum(sum(Swkf))*df*dk));

    % Amplitude for wave(km, fn)
    Umn = sqrt(2*Sukf*dk*df);
    Vmn = sqrt(2*Svkf*dk*df);
    Wmn = sqrt(2*Swkf*dk*df);

    % Divergence-free condition
    kxmn = -kvec .* (Vmn ./ Umn .* cos(PSI) + Wmn ./ Umn .* sin(PSI));
    kxmn(:,1) = 0;

    % Phase for wave(fn)
    Phmn = kxmn*X(i)+kvec.*cos(PSI)*Y(i)+kvec.*sin(PSI)*Z(i)+PHI;

    % FFT matrix of u time series
    Uf(1:N,:) = (2*N-1)*complex(1/2*Umn.*cos(Phmn), 1/2*Umn.*sin(Phmn))';
    Uf(2*N-1:-1:N+1,:) = ...
        (2*N-1)*complex(1/2*Umn(:,2:end).*cos(-Phmn(:,2:end)),...
        1/2*Umn(:,2:end).*sin(-Phmn(:,2:end)))';

    % FFT matrix of v time series
    Vf(1:N,:) = (2*N-1)*complex(1/2*Vmn.*cos(Phmn), 1/2*Vmn.*sin(Phmn))';
    Vf(2*N-1:-1:N+1,:) = ...
        (2*N-1)*complex(1/2*Vmn(:,2:end).*cos(-Phmn(:,2:end)), ...
        1/2*Vmn(:,2:end).*sin(-Phmn(:,2:end)))';

    % FFT matrix of w time series
    Wf(1:N,:) = (2*N-1)*complex(1/2*Wmn.*cos(Phmn), 1/2*Wmn.*sin(Phmn))';
    Wf(2*N-1:-1:N+1,:) = ...
        (2*N-1)*complex(1/2*Wmn(:,2:end).*cos(-Phmn(:,2:end)), ...
        1/2*Wmn(:,2:end).*sin(-Phmn(:,2:end)))';

    % Inverse FFT to calculate the velocity time series
    u(:,i) = sum(ifft(Uf), 2);
    v(:,i) = sum(ifft(Vf), 2);
    w(:,i) = sum(ifft(Wf), 2);

    fprintf('Pt. %3d completed\n\n', i);
end

%% Save data to HDF5 database
h5write(hdf5File,'/GRID',GRID');
h5write(hdf5File,'/TIME',tvec);
h5write(hdf5File,'/UMEAN',Uav);
h5write(hdf5File,'/U',u);
h5write(hdf5File,'/V',v);
h5write(hdf5File,'/W',w);

%% Post-processing scripts
post_psd
post_coh

