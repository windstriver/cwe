%% Run parameters.m to set primary variables
run('parameters');

%% Generate random phase matrices
PSI = 2*pi*rand(M,N);
PHI = 2*pi*rand(M,N);

%% Start generating turbulence

% Velocity time series
% row: time step index
% column: points index
u = zeros(nt,nd);
v = zeros(nt,nd);
w = zeros(nt,nd);

% Initialization
% Wavenumber-Freq Spectrum matrix at a single point
Sukf = zeros(M,N);    % Suu(km, fn)
Svkf = zeros(M,N);    % Svv(km, fn)
Swkf = zeros(M,N);    % Sww(km, fn)
% Amplitude for wave(km, fn)
Umn = zeros(M,N);
Vmn = zeros(M,N);
Wmn = zeros(M,N);
% kxmn to maintain divergence-free condition
kxmn = zeros(M,N);
% Phase for wave(fn)
Phmn = zeros(M,N);

%% Open the HDF5 database
fileattrib(hdf5File,'+w');
plist = 'H5P_DEFAULT';
fid = H5F.open(hdf5File,'H5F_ACC_RDWR',plist);
% write dataset /GRID
dset_id = H5D.open(fid,'/GRID');
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,GRID');
% write dataset /TIME
dset_id = H5D.open(fid,'/TIME');
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,tvec);
% write dataset /UMEAN
dset_id = H5D.open(fid,'/UMEAN');
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,Uav);

% read the FLAG dataset to determine whether the simulation status
dset_id = H5D.open(fid,'/FLAG');
% Flag to store the simulation status
FLAG = H5D.read(dset_id);

% parpool to specify workers
c = parcluster('local');
c.NumWorkers = 35;
parpool(c, c.NumWorkers);

% start and end index of the points to be simulated
ptStart = 1;
ptEnd = 5000;

parfor i = ptStart:ptEnd
    % i: points index
    % FFT matrix for velocity time series u, v, w
    % row: frequency index
    % 0, df, 2*df, ... , (N-1)*df
    % N*df, (N+1)*df, ... , (2*N-2)*df
    % column: wavenumber index
    Uf = zeros(2*N-1, M);
    Vf = zeros(2*N-1, M);
    Wf = zeros(2*N-1, M);
    % Wave-freq. Spectrum matrix for point i
    Sukf = 2/pi * ones(M,1) * (Su0(:,i)'.*(Cxyz(1)*fvec/U_ref)) ./ ...
        ((Cxyz(1)*fvec/U_ref).^2 + kvec.^2);
    Svkf = 2/pi * ones(M,1) * (Sv0(:,i)'.*(Cxyz(2)*fvec/U_ref)) ./ ...
        ((Cxyz(2)*fvec/U_ref).^2 + kvec.^2);
    Swkf = 2/pi * ones(M,1) * (Sw0(:,i)'.*(Cxyz(3)*fvec/U_ref)) ./ ...
        ((Cxyz(3)*fvec/U_ref).^2 + kvec.^2);

    % Compare std. calculated by integration of wave-freq. spectrum with
    % theory
    fprintf('Pt. %3d (%6.3f, %6.3f, %6.3f):\n', i, X(i), Y(i), Z(i));

    urmsAuto = Iu(i)*Uav(i);
    urmsWave = sqrt(sum(sum(Sukf))*df*dk);
    fprintf('std(u) = %6.3f by auto-spectra\n', urmsAuto);
    fprintf('std(u) = %6.3f (%4.2f) by wave-freq spectra\n',...
            urmsWave, 100*abs(urmsWave-urmsAuto)/urmsAuto);

    vrmsAuto = Iv(i)*Uav(i);
    vrmsWave = sqrt(sum(sum(Svkf))*df*dk);
    fprintf('std(v) = %6.3f by auto-spectra\n', vrmsAuto);
    fprintf('std(v) = %6.3f (%4.2f) by wave-freq spectra\n',...
            vrmsWave, 100*abs(vrmsWave-vrmsAuto)/vrmsAuto);

    wrmsAuto = Iw(i)*Uav(i);
    wrmsWave = sqrt(sum(sum(Swkf))*df*dk);
    fprintf('std(w) = %6.3f by auto-spectra\n', wrmsAuto);
    fprintf('std(w) = %6.3f (%4.2f) by wave-freq spectra\n',...
            wrmsWave, 100*abs(wrmsWave-wrmsAuto)/wrmsAuto);

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

    fprintf('Pt. %6d(%3d)completed\n\n', i, nd);

    FLAG(i) = 1;
end

%% Save data to HDF5 database
% HDF5 block settings
blockStart = [0 ptStart-1];
h5_start = fliplr(blockStart);
block = [nt ptEnd-ptStart+1];
h5_block = fliplr(block);
% write dataset /U
dset_id = H5D.open(fid,'/U');
mem_space_id = H5S.create_simple(2,h5_block,[]);
file_space_id = H5D.get_space(dset_id);
H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',...
    h5_start,[],[],h5_block);
H5D.write(dset_id,'H5ML_DEFAULT',...
    mem_space_id,file_space_id,plist,u(:,ptStart:ptEnd));
% write dataset /V
dset_id = H5D.open(fid,'/V');
mem_space_id = H5S.create_simple(2,h5_block,[]);
file_space_id = H5D.get_space(dset_id);
H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',...
    h5_start,[],[],h5_block);
H5D.write(dset_id,'H5ML_DEFAULT',...
    mem_space_id,file_space_id,plist,v(:,ptStart:ptEnd));
% write dataset /W
dset_id = H5D.open(fid,'/W');
mem_space_id = H5S.create_simple(2,h5_block,[]);
file_space_id = H5D.get_space(dset_id);
H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',...
    h5_start,[],[],h5_block);
H5D.write(dset_id,'H5ML_DEFAULT',...
    mem_space_id,file_space_id,plist,w(:,ptStart:ptEnd));

% write dataset /FLAG
dset_id = H5D.open(fid,'/FLAG');
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,FLAG);

% close dataset and hdf5 file
H5D.close(dset_id);
H5F.close(fid);

fprintf('\n\nSimulation finished for pt. %6d to %6d \n', ptStart, ptEnd);

