function scf = spatialCorrelationFactor(d,C,U,f,Sui,Suj)
%SPATIALCORRELATIONFACTOR Calculate the spatial correlation factor of two
% random process ui(t)(at node i) and uj(t)(at node j)
% Input:
% --d: distance between node i and node j
% --C: decay factor in coherency function
% --U: mean velocity
% --f: frequency vector
% --Sui: auto-spectra of random process ui(t)
% --Suj: auto-spectra of random process uj(t)
% Output:
% --scf: spatial correlation factor
coh = exp(-C*f'*d/U);
scf = trapz(sqrt(Sui.*Suj).*coh) / ...
    (sqrt(trapz(Sui))*sqrt(trapz(Suj)));
end
