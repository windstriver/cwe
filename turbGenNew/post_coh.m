%% Comparison of coherency function with that of Davenport
[Coh, freqCoh] = mscohere(u(:,2),u(:,4),hann(512),256,512,1/dt);
fvecRe = fvec*(Z(4)-Z(2))/((Uav(2)+Uav(4))/2);
figure;
plot(freqCoh,Coh,'b-',...
    fvec,exp(-Cxyz(1)*fvecRe),'r--')
xlabel('$ f(\mathrm{Hz}) $', 'Interpreter', 'latex')
ylabel('Coherency')
legend('Simulation', 'Target')
xlim([0,25])
title('Coherency function');
% saveas(gcf, 'coh', 'epsc')

%% Spatial correlation factor
% Target
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
ylabel('Spatial Correlation Factor');
legend('Target','Simulation');
title('Spatial Correlation Factor');
