h5file = 'inflowTurb.h5';
h5file2 = 'inflowTurbInt.h5';
% Scaling for turbulence intensity
% sigma_scale = 11.1438 * 0.116;
sigma_scale = 1;

delta_t = h5readatt(h5file, '/', 'delta_x1');
delta_y = h5readatt(h5file, '/', 'delta_x2');
delta_z = h5readatt(h5file, '/', 'delta_x3');

M1 = h5readatt(h5file, '/', 'M1');
M2 = h5readatt(h5file, '/', 'M2');
M3 = h5readatt(h5file, '/', 'M3');

Lt0 = M1 * delta_t;
Ly0 = M2 * delta_y;
Lz0 = M3 * delta_z;

[y_sr, z_sr] = ndgrid((0:M2-1)*delta_y, (0:M3-1)*delta_z);

% CFD inlet face centres
face_centres = readmatrix('faceCentresInlet.csv');
% reorder coordinates
y_grid = -face_centres(:,4);
z_grid = face_centres(:,3);

% simulation time
deltaT = 5e-4;
endT = Lt0;
t_vec = (0:deltaT:endT)';

% Modification using periodic property
y_grid_mod = mod(y_grid, Ly0);
z_grid_mod = mod(z_grid, Lz0);
t_vec_mod = mod(t_vec, Lt0);

%% Interpolation
nd = length(y_grid);
nt = length(t_vec);

if exist(h5file2)
    disp('HDF5 file already exists');
else
    fid = H5F.create(h5file2);
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
    H5D.create(fid,'V1',type_id,space_id,dcpl);
    H5D.create(fid,'V3',type_id,space_id,dcpl);
    H5D.create(fid,'W',type_id,space_id,dcpl);

    H5S.close(space_id);
    H5T.close(type_id);
    H5F.close(fid);
end

%% Mean velocity
H = 0.5;
U_ref = 11.1438;
delta = 2*H;
u_mean = U_ref .* (z_grid/H).^(1/4);
u_mean(z_grid>delta) = U_ref * (delta/H).^(1/4);
h5write(h5file2, '/UMEAN', u_mean);

%% Save attributes to HDF5
h5writeatt(h5file, '/', 'deltaT', deltaT);
h5writeatt(h5file, '/', 'endT', endT);

h5create(h5file2, '/GRID', size(face_centres(:,3:4)));
h5write(h5file2, '/GRID', face_centres(:,3:4));

for i = 1:nt
    if t_vec(i) == Lt0
        t_vec(i) = t_vec(i) - 1e-3*delta_t;
    end
    
    i_l = find(~(t_vec(i) - (0:M1-1)*delta_t >= 0), 1) - 1;
    
    % interpolation of u
    u_sr = h5read(h5file, '/u', [i_l 1 1], [2 Inf Inf]);
    u_sr_int1 = griddedInterpolant(y_sr, z_sr, squeeze(u_sr(1,:,:)));
    u_sr_int2 = griddedInterpolant(y_sr, z_sr, squeeze(u_sr(2,:,:)));
    u_int1 = u_sr_int1(y_grid_mod, z_grid_mod);
    u_int2 = u_sr_int2(y_grid_mod, z_grid_mod);
    u_int = (t_vec_mod(i)-(i_l-1)*delta_t)/delta_t * u_int1 + ...
        (i_l*delta_t - t_vec_mod(i))/delta_t * u_int2;
    
    h5write(h5file2, '/U', sigma_scale * u_int', [i 1], [1 nd]);
    
    % interpolation of v and save to w (reorder)
    u_sr = h5read(h5file, '/v', [i_l 1 1], [2 Inf Inf]);
    u_sr_int1 = griddedInterpolant(y_sr, z_sr, squeeze(u_sr(1,:,:)));
    u_sr_int2 = griddedInterpolant(y_sr, z_sr, squeeze(u_sr(2,:,:)));
    u_int1 = u_sr_int1(y_grid_mod, z_grid_mod);
    u_int2 = u_sr_int2(y_grid_mod, z_grid_mod);
    u_int = (t_vec_mod(i)-(i_l-1)*delta_t)/delta_t * u_int1 + ...
        (i_l*delta_t - t_vec_mod(i))/delta_t * u_int2;
    
    h5write(h5file2, '/W', -sigma_scale * u_int', [i 1], [1 nd]);
    
    % interpolation of w
    % first sample
    u_sr = h5read(h5file, '/w1', [i_l 1 1], [2 Inf Inf]);
    u_sr_int1 = griddedInterpolant(y_sr, z_sr, squeeze(u_sr(1,:,:)));
    u_sr_int2 = griddedInterpolant(y_sr, z_sr, squeeze(u_sr(2,:,:)));
    u_int1 = u_sr_int1(y_grid_mod, z_grid_mod);
    u_int2 = u_sr_int2(y_grid_mod, z_grid_mod);
    u_int = (t_vec_mod(i)-(i_l-1)*delta_t)/delta_t * u_int1 + ...
        (i_l*delta_t - t_vec_mod(i))/delta_t * u_int2;
    
    h5write(h5file2, '/V1', sigma_scale * u_int', [i 1], [1 nd]);

    % third sample
    u_sr = h5read(h5file, '/w3', [i_l 1 1], [2 Inf Inf]);
    u_sr_int1 = griddedInterpolant(y_sr, z_sr, squeeze(u_sr(1,:,:)));
    u_sr_int2 = griddedInterpolant(y_sr, z_sr, squeeze(u_sr(2,:,:)));
    u_int1 = u_sr_int1(y_grid_mod, z_grid_mod);
    u_int2 = u_sr_int2(y_grid_mod, z_grid_mod);
    u_int = (t_vec_mod(i)-(i_l-1)*delta_t)/delta_t * u_int1 + ...
        (i_l*delta_t - t_vec_mod(i))/delta_t * u_int2;
    
    h5write(h5file2, '/V3', sigma_scale * u_int', [i 1], [1 nd]);


    fprintf('Time index %u finished.\n', i);
end

