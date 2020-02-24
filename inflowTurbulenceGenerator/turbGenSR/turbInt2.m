h5file = '../inflowTurb8.h5';
h5file2 = 'inflowTurbInt8.h5';

% Scaling for turbulence intensity
sigma_scale = 11.1438 * 0.116;

delta_t = h5readatt(h5file, '/', 'delta_x1');
delta_y = h5readatt(h5file, '/', 'delta_x2');
delta_z = h5readatt(h5file, '/', 'delta_x3');

M1 = h5readatt(h5file, '/', 'M1');
M2 = h5readatt(h5file, '/', 'M2');
M3 = h5readatt(h5file, '/', 'M3');

Lt0 = M1 * delta_t;
Ly0 = M2 * delta_y;
Lz0 = M3 * delta_z;

% [y_sr, z_sr] = ndgrid((0:M2-1)*delta_y, (0:M3-1)*delta_z);
t_sr = (0:M1-1)*delta_t;
y_sr = (0:M2-1)*delta_y;
z_sr = (0:M3-1)*delta_z;

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

%% Mean velocity
H = 0.5;
U_ref = 11.1438;
delta = 2*H;
u_mean = U_ref .* (z_grid/H).^(1/4);
u_mean(z_grid>delta) = U_ref * (delta/H).^(1/4);

%% Interpolation
nd = length(y_grid);
nt = length(t_vec);

if exist(h5file2, 'file')
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
    
    h5_dims = [2 nd];
    h5_maxdims = h5_dims;
    space_id = H5S.create_simple(2,h5_dims,h5_maxdims);
    dcpl = 'H5P_DEFAULT';
    H5D.create(fid,'GRID',type_id,space_id,dcpl);

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

%% Save to HDF5
h5writeatt(h5file, '/', 'deltaT', deltaT);
h5writeatt(h5file, '/', 'endT', endT);
h5write(h5file2, '/GRID', face_centres(:,3:4));
h5write(h5file2, '/UMEAN', u_mean);

for i = 1:nd
    % find location closest to ith pt
    [~, index_y] = min(abs(y_grid_mod(i) - y_sr));
    [~, index_z] = min(abs(z_grid_mod(i) - z_sr));
    
    % interpolation of u
    u_sr = h5read(h5file, '/u', [1 index_y index_z], [M1 1 1]);
    u_sr_int = griddedInterpolant(t_sr', u_sr);
    u_int = u_sr_int(t_vec_mod);
    
    h5write(h5file2, '/U', sigma_scale * u_int, [1 i], [nt 1]);
    
    % interpolation of v and save to W (reorder)
    u_sr = h5read(h5file, '/v', [1 index_y index_z], [M1 1 1]);
    u_sr_int = griddedInterpolant(t_sr', u_sr);
    u_int = u_sr_int(t_vec_mod);
    
    h5write(h5file2, '/W', -sigma_scale * u_int, [1 i], [nt 1]);
    
    % interpolation of w and save to V (reorder)
    u_sr = h5read(h5file, '/w1', [1 index_y index_z], [M1 1 1]);
    u_sr_int = griddedInterpolant(t_sr', u_sr);
    u_int = u_sr_int(t_vec_mod);
    
    h5write(h5file2, '/V1', sigma_scale * u_int, [1 i], [nt 1]);
    
    % interpolation of w and save to V (reorder)
    u_sr = h5read(h5file, '/w3', [1 index_y index_z], [M1 1 1]);
    u_sr_int = griddedInterpolant(t_sr', u_sr);
    u_int = u_sr_int(t_vec_mod);
    
    h5write(h5file2, '/V3', sigma_scale * u_int, [1 i], [nt 1]);
    
    fprintf('Position index %u finished.\n', i);
end
