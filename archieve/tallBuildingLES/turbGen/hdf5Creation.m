%% Create HDF5 database to store simulated data
% HDF5 library uses C-style ordering for multidimensional arrays
% MATLAB uses FORTRAN-style ordering
% Dataset size in HDF5 database for use in OpenFOAM
% GRID: [nd 4]    PtNum, x, y, z
% TIME: [nt]
% UMEAN: [nd]
% U,V,W: [nd nt]
fid = H5F.create(hdf5File);
type_id = H5T.copy('H5T_NATIVE_DOUBLE');
order = H5ML.get_constant_value('H5T_ORDER_BE');
H5T.set_order(type_id,order);

% dataset /GRID
h5_dims = [nd 4];
h5_maxdims = h5_dims;
space_id = H5S.create_simple(2,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'GRID',type_id,space_id,dcpl);

% dataset /TIME
h5_dims = nt;
h5_maxdims = h5_dims;
space_id = H5S.create_simple(1,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'TIME',type_id,space_id,dcpl);

% dataset /UMEAN
h5_dims = nd;
h5_maxdims = h5_dims;
space_id = H5S.create_simple(1,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'UMEAN',type_id,space_id,dcpl);

% dataset /FLAG
h5_dims = nd;
h5_maxdims = h5_dims;
space_id = H5S.create_simple(1,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'FLAG',type_id,space_id,dcpl);

% dataset /U
h5_dims = [nd nt];
h5_maxdims = h5_dims;
space_id = H5S.create_simple(2,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'U',type_id,space_id,dcpl);

% dataset /V
h5_dims = [nd nt];
h5_maxdims = h5_dims;
space_id = H5S.create_simple(2,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'V',type_id,space_id,dcpl);

% dataset /W
h5_dims = [nd nt];
h5_maxdims = h5_dims;
space_id = H5S.create_simple(2,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'W',type_id,space_id,dcpl);

H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);
H5F.close(fid);

%h5disp(hdf5File);