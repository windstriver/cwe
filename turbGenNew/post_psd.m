%% Along-wind speed spectra
[Suu,freq] = pwelch(u(:,10),hann(2048),1024,2048,1/dt,'onesided');
figure;
loglog(freq,Suu,'b-',fvec,Su0(:,10),'r--');
xlim([1,100]);
title('Along-wind speed spectra');
xlabel('$ f, \mathrm{Hz} $', 'Interpreter', 'latex');
ylabel('$ S_{uu}, \mathrm{m^2 s^2 / Hz} $', 'Interpreter', 'latex');
legend('Simulation','Target');
title('Along-wind spectra');
saveas(gcf, 'Suu', 'epsc')

%% Cross-wind speed spectra
[Svv,freq] = pwelch(v(:,10),hann(2048),1024,2048,1/dt,'onesided');
figure;
loglog(freq,Svv,'b-',fvec,Sv0(:,10),'r--');
xlim([1,100]);
title('Cross-wind speed spectra');
xlabel('$ f, \mathrm{Hz} $', 'Interpreter', 'latex');
ylabel('$ S_{vv}, \mathrm{m^2 s^2 / Hz} $', 'Interpreter', 'latex');
legend('Simulation','Target');
title('Cross-wind spectra');
saveas(gcf, 'Svv', 'epsc')

%% Vertical wind speed spectra
[Sww,freq] = pwelch(w(:,10),hann(2048),1024,2048,1/dt,'onesided');
figure;
loglog(freq,Sww,'b-',fvec,Sw0(:,10),'r--');
xlim([1,100]);
title('Vertical wind speed spectra');
xlabel('$ f, \mathrm{Hz} $', 'Interpreter', 'latex');
ylabel('$ S_{ww}, \mathrm{m^2 s^2 / Hz} $', 'Interpreter', 'latex');
legend('Simulation','Target');
title('Vertical wind spectra');
saveas(gcf, 'Sww', 'epsc')

