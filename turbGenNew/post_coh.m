%% Read turbulent velocity from HDF5 database
u = h5read(hdf5File, '/U');
v = h5read(hdf5File, '/V');
w = h5read(hdf5File, '/W');

%% Comparison of coherency function with that of Davenport
pt1 = find(Z==0.1);    % point z = 0.1 index
pt2 = find(Z==0.3);    % point z = 0.3 index
% Coherency of u velocity between pt1 and pt2
[msCohu, freqCohu] = mscohere(u(:,pt1), u(:,pt2), hann(512), 256, 512, 1/dt);
[msCohv, freqCohv] = mscohere(v(:,pt1), v(:,pt2), hann(512), 256, 512, 1/dt);
[msCohw, freqCohw] = mscohere(w(:,pt1), w(:,pt2), hann(512), 256, 512, 1/dt);

figure;
plot(freqCohu, msCohu, 'b--', ...
    freqCohv, msCohv, 'g:', ...
    freqCohw, msCohw, 'm-.', ...
    fvec, exp(-Cxyz(1)*fvec*(Z(pt2)-Z(pt1))/((Uav(pt1)+Uav(pt2))/2)),'r-')
xlabel('$ f(\mathrm{Hz}) $', 'Interpreter', 'latex');
ylabel('Coherency');
leg1 = legend('$Coh(u_1, u_2)$', '$Coh(v_1, v_2)$', '$Coh(w_1, w_2)$', 'Target');
set(leg1, 'Interpreter', 'latex');
xlim([0,25]);
title('Coherency of velocity between two points $(z_1=0.1\,\mathrm{m}, z_2=0.3\,\mathrm{m})$', ...
    'Interpreter', 'latex');
saveas(gcf, 'homo_turb_coh_uvw', 'svg');

%% Spatial correlation factor
dzvec = zeros(1,length(Z));
scfTarget = zeros(1,length(Z));
scfSim = zeros(1,length(Z));
scfTarget(1) = spatialCorrelationFactor(0,Cxyz(1),...
        Uav(1),fvec,Su0(:,1),Su0(:,1));
R = corrcoef(u(:,1),u(:,1));
scfSim(1) = R(1,2);
for i = 2:length(Z)
    dzvec(i) = Z(i) - Z(1);
    scfTarget(i) = spatialCorrelationFactor(dzvec(i),Cxyz(1),...
        (Uav(1)+Uav(i))/2,fvec,Su0(:,1),Su0(:,i));
    R = corrcoef(u(:,1),u(:,i));
    scfSim(i) = R(1,2);
end
figure;
plot(dzvec,scfTarget,'r--',dzvec,scfSim,'b-');
xlabel('Vertical distance (m)');
ylabel('$Sc_{i,j}$', 'Interpreter', 'latex');
legend('Target','Simulation');
title('Correlation coefficient of $u_i$ and $u_j$', ...
    'Interpreter', 'latex');
saveas(gcf, 'homo_turb_scf', 'svg')
