%% Mass flux correction
h5fileInt = 'inflowTurbInt.h5';
h5fileMFC = 'inflowTurbIntMFC.h5';

u_mean = h5read(h5fileInt, '/UMEAN');

nd = length(u_mean);

deltaT = 5e-4;
endT = 4.8+deltaT;

t_vec = (0:deltaT:endT)';
nt = length(t_vec);

face_area = readmatrix('faceAreaVectorsInlet.csv');
Sf = -face_area(:,1)';
S = sum(Sf);
Ub0 = sum(u_mean'.*Sf) / S;

% create HDF5 with MFC
if exist(h5fileMFC)
    disp('HDF5 file already exists')
else
    fid = H5F.create(h5fileMFC);
    type_id = H5T.copy('H5T_NATIVE_DOUBLE');
    order = H5ML.get_constant_value('H5T_ORDER_BE');
    H5T.set_order(type_id,order);

    h5_dims = nd;
    h5_maxdims = h5_dims;
    space_id = H5S.create_simple(1,h5_dims,h5_maxdims);
    dcpl = 'H5P_DEFAULT';
    dset_id = H5D.create(fid,'UMEAN',type_id,space_id,dcpl);

    h5_dims = [nd nt];
    h5_maxdims = h5_dims;
    space_id = H5S.create_simple(2,h5_dims,h5_maxdims);
    dcpl = 'H5P_DEFAULT';
    H5D.create(fid,'U',type_id,space_id,dcpl);
    H5D.create(fid,'V',type_id,space_id,dcpl);
    H5D.create(fid,'W',type_id,space_id,dcpl);

    H5S.close(space_id);
    H5T.close(type_id);
    H5F.close(fid);
end

h5write(h5fileMFC, '/UMEAN', u_mean);

%% MFC
Ubt = zeros(nt,1);
for i = 1:nt
    u = h5read(h5fileInt, '/U', [i 1], [1 Inf]);
    v = h5read(h5fileInt, '/V', [i 1], [1 Inf]);
    w = h5read(h5fileInt, '/W', [i 1], [1 Inf]);

    Ubt(i) = sum(u.*Sf) / sum(Sf) + Ub0;
    u_mfc = Ub0 / Ubt(i) * u;
    v_mfc = Ub0 / Ubt(i) * v;
    w_mfc = Ub0 / Ubt(i) * w;
    
    h5write(h5fileMFC, '/U', u_mfc, [i 1], [1 nd]);
    h5write(h5fileMFC, '/V', v_mfc, [i 1], [1 nd]);
    h5write(h5fileMFC, '/W', w_mfc, [i 1], [1 nd]);
end

mfc_ratio = Ub0 ./ Ubt;
save('mfc_ratio.txt', 'mfc_ratio', '-ascii');

