%% Run parameters to set primary variables
run('parameters.m');

%% Read simulated velocity from HDF5 database
u = h5read(hdf5File,'/U');

%% Spatial correlation factor
% assume only mean velocity follows power law
% turbulence intensity and length scale are constant
% stdRatio to measure the degree of inhomogenity
% stdRatio = std(u_i) / std(u_ref) - 1
% where u_i is the velocity time history at pt i
% u_ref is the velocity time history at reference point
stdRatio = (Z/Z(4)).^alphau - 1;
scfTarget = zeros(1,length(Z));
scfSim = zeros(1,length(Z));

for i = 1:length(Z)
    scfTarget(i) = spatialCorrelationFactor(abs(Z(i)-Z(4)),Cxyz(1),...
        (Uav(i)+Uav(4))/2,fvec,Su0(:,4),Su0(:,4));
    R = corrcoef(u(:,4),u(:,i));
    scfSim(i) = R(1,2);
end

%% Plot spatial correlation factor versus std. ratio
figure;
plot(stdRatio*100,scfTarget,'r--',stdRatio*100,scfSim,'b-');
xlabel('$(\sigma_{ui}-\sigma_{uref})/\sigma_{uref}$',...
       'Interpreter','latex');
ylabel('Spatial correlation factor');
title('Spatial Correlation Factor');
legend('Target','Simulation');
saveas(gcf,'scf_inhomo','svg');
