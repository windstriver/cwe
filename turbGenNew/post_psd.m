%% Run parameters to set primary variables
run('parameters.m');

%% Read simulated velocity from HDF5 database
hdf5File = 'homo_turb_validation.h5';
u = h5read(hdf5File,'/U');
v = h5read(hdf5File,'/V');
w = h5read(hdf5File,'/W');

%% PSD
pt1 = find(Z==0.1);    % point z = 0.1 index
pt2 = find(Z==0.5);    % point z = 0.5 index
% Along-wind speed spectra
[Suu1,freq1] = pwelch(u(:,pt1),hann(2048),1024,2048,1/dt,'onesided');
[Suu2,freq2] = pwelch(u(:,pt2),hann(2048),1024,2048,1/dt,'onesided');
figure;
loglog(freq1,Suu1,'b-',fvec,Su0(:,pt1),'r--', ...
    freq2,Suu2,'g-',fvec,Su0(:,pt2),'r-.');
xlim([1,100]);
xlabel('$ f, \mathrm{Hz} $', 'Interpreter', 'latex');
ylabel('$ S_{uu}, \mathrm{m^2 s^2 / Hz} $', 'Interpreter', 'latex');
leg1 = legend('Simulation $z=0.1\,\mathrm{m}$','Target $z=0.1\,\mathrm{m}$', ...
    'Simulation $z=0.5\,\mathrm{m}$','Target $z=0.5\,\mathrm{m}$');
set(leg1, 'Interpreter', 'latex');
title('Along-wind spectrum');
saveas(gcf, 'homo_psd_uu', 'svg')

% Cross-wind speed spectra
[Svv1,freq1] = pwelch(v(:,pt1),hann(2048),1024,2048,1/dt,'onesided');
[Svv2,freq2] = pwelch(v(:,pt2),hann(2048),1024,2048,1/dt,'onesided');
figure;
loglog(freq1,Svv1,'b-',fvec,Sv0(:,pt1),'r--', ...
    freq2,Svv2,'g-',fvec,Sv0(:,pt2),'r-.');
xlim([1,100]);
title('Cross-wind speed spectrum');
xlabel('$ f, \mathrm{Hz} $', 'Interpreter', 'latex');
ylabel('$ S_{vv}, \mathrm{m^2 s^2 / Hz} $', 'Interpreter', 'latex');
leg1 = legend('Simulation $z=0.1\,\mathrm{m}$','Target $z=0.1\,\mathrm{m}$', ...
    'Simulation $z=0.5\,\mathrm{m}$','Target $z=0.5\,\mathrm{m}$');
set(leg1, 'Interpreter', 'latex');
title('Cross-wind spectrum');
saveas(gcf, 'homo_psd_vv', 'svg')

% Vertical wind speed spectra
[Sww1,freq1] = pwelch(w(:,pt1),hann(2048),1024,2048,1/dt,'onesided');
[Sww2,freq2] = pwelch(w(:,pt2),hann(2048),1024,2048,1/dt,'onesided');
figure;
loglog(freq1,Sww1,'b-',fvec,Sw0(:,pt1),'r--', ...
    freq1,Sww2,'g-',fvec,Sw0(:,pt2),'r-.');
xlim([1,100]);
xlabel('$ f, \mathrm{Hz} $', 'Interpreter', 'latex');
ylabel('$ S_{ww}, \mathrm{m^2 s^2 / Hz} $', 'Interpreter', 'latex');
leg1 = legend('Simulation $z=0.1\,\mathrm{m}$','Target $z=0.1\,\mathrm{m}$', ...
    'Simulation $z=0.5\,\mathrm{m}$','Target $z=0.5\,\mathrm{m}$');
set(leg1, 'Interpreter', 'latex');
title('Vertical wind speed spectrum');
saveas(gcf, 'homo_psd_ww', 'svg')
