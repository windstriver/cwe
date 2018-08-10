%% Run parameters to set primary variables
% run('parameters.m');

%% Read simulated velocity from HDF5 database
u = h5read(hdf5File,'/U');
v = h5read(hdf5File,'/V');
w = h5read(hdf5File,'/W');
[~, pt0] = min(abs(Uav-10));

%% Mean velocity
figure;
plot(Uav, Z);
hold on;
plot([10, 10], [0, Z(pt0)], 'r--');
plot([5, 10], [Z(pt0) Z(pt0)], 'r--');
text(6, 0.9, '$U_{av}=U_{avref}\left(\frac{z}{z_{ref}}\right)^\alpha$',...
     'Interpreter', 'latex');
text(6, 0.8, '$U_{avref} = 10\,\mathrm{m/s},\, z_{ref}=0.364\,\mathrm{m},\,\alpha=0.326$',...
     'Interpreter', 'latex');
hold off;
xlabel('$U_{av} (\mathrm{m/s})$', 'Interpreter', 'latex');
ylabel('$z (\mathrm{m})$', 'Interpreter', 'latex');
title('Mean velocity profile')
saveas(gcf,'inhomo_Uav','svg');

%% Turbulent intensity
figure;
plot(100*Iu, Z, 'r-', 100*Iv, Z, 'b--', 100*Iw, Z, 'g-.');
text(22, 0.9, '$I_{j}=I_{refj}\left(\frac{z}{z_{ref}}\right)^{-d_j}$',...
     'Interpreter', 'latex');
text(22, 0.8, '$I_{refj} = 0.208, 0.182, 0.152$',...
     'Interpreter', 'latex');
text(22, 0.7, '$d_j = 0.191, 0.123, 0.005$',...
     'Interpreter', 'latex');
text(22, 0.6, 'in the $u, v, w$ directions, respectively',...
     'Interpreter', 'latex');
xlabel('$I(\%)$', 'Interpreter', 'latex');
ylabel('$z (\mathrm{m})$', 'Interpreter', 'latex');
leg1 = legend('$I_u$', '$I_v$', '$I_w$');
set(leg1, 'Interpreter', 'latex');
title('Turbulence intensity profile');
saveas(gcf,'inhomo_TI','svg');

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
saveas(gcf, 'inhomo_psd_uu', 'svg')

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
saveas(gcf, 'inhomo_psd_vv', 'svg')

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
saveas(gcf, 'inhomo_psd_ww', 'svg')

%% Spatial correlation factor
scfTarget = zeros(1,length(Z));
scfSim = zeros(1,length(Z));
for i = 1:length(Z)
    scfTarget(i) = spatialCorrelationFactor(abs(Z(i)-Z(pt0)),Cxyz(1),...
        (Uav(i)+Uav(pt0))/2,fvec,Su0(:,pt0),Su0(:,pt0));
    R = corrcoef(u(:,pt0),u(:,i));
    scfSim(i) = R(1,2);
end

% Plot spatial correlation factor versus std. ratio
figure;
plot(scfTarget, Z, 'r--',scfSim, Z, 'b-');
xlabel('$Sc_{i,0}$',...
       'Interpreter','latex');
ylabel('$z (\mathrm{m})$', 'Interpreter', 'latex');
title('Correlation coefficient of $u_i$ and $u_0$', ...
    'Interpreter', 'latex');
legend('Target','Simulation');
saveas(gcf,'inhomo_scf','svg');
