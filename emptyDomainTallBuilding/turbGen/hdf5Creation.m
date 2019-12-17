%% Create HDF5 database to store simulated data
% HDF5 library uses C-style h5.ordering for multidimensional arrays
% MATLAB uses FORTRAN-style h5.ordering
% Dataset size in HDF5 database for use in OpenFOAM
% GRID: [nd 4]    PtNum, x, y, z
% TIME: [nt]
% UMEAN: [nd]
% U,V,W: [nd nt]
h5.fid = H5F.create(hdf5File);
h5.type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5.order = H5ML.get_constant_value('H5T_ORDER_BE');
H5T.set_order(h5.type_id,h5.order);

% dataset /GRID
h5.dims = [nd 4];
h5.maxdims = h5.dims;
h5.space_id = H5S.create_simple(2,h5.dims,h5.maxdims);
h5.dcpl = 'H5P_DEFAULT';
h5.dset_id = H5D.create(h5.fid,'GRID',h5.type_id,h5.space_id,h5.dcpl);

% dataset /TIME
h5.dims = nt;
h5.maxdims = h5.dims;
h5.space_id = H5S.create_simple(1,h5.dims,h5.maxdims);
h5.dcpl = 'H5P_DEFAULT';
h5.dset_id = H5D.create(h5.fid,'TIME',h5.type_id,h5.space_id,h5.dcpl);

% dataset /UMEAN
h5.dims = nd;
h5.maxdims = h5.dims;
h5.space_id = H5S.create_simple(1,h5.dims,h5.maxdims);
h5.dcpl = 'H5P_DEFAULT';
h5.dset_id = H5D.create(h5.fid,'UMEAN',h5.type_id,h5.space_id,h5.dcpl);

% dataset /FLAG
h5.dims = nd;
h5.maxdims = h5.dims;
h5.space_id = H5S.create_simple(1,h5.dims,h5.maxdims);
h5.dcpl = 'H5P_DEFAULT';
h5.dset_id = H5D.create(h5.fid,'FLAG',h5.type_id,h5.space_id,h5.dcpl);

% dataset /U
h5.dims = [nd nt];
h5.maxdims = h5.dims;
h5.space_id = H5S.create_simple(2,h5.dims,h5.maxdims);
h5.dcpl = 'H5P_DEFAULT';
h5.dset_id = H5D.create(h5.fid,'U',h5.type_id,h5.space_id,h5.dcpl);

% dataset /V
h5.dims = [nd nt];
h5.maxdims = h5.dims;
h5.space_id = H5S.create_simple(2,h5.dims,h5.maxdims);
h5.dcpl = 'H5P_DEFAULT';
h5.dset_id = H5D.create(h5.fid,'V',h5.type_id,h5.space_id,h5.dcpl);

% dataset /W
h5.dims = [nd nt];
h5.maxdims = h5.dims;
h5.space_id = H5S.create_simple(2,h5.dims,h5.maxdims);
h5.dcpl = 'H5P_DEFAULT';
h5.dset_id = H5D.create(h5.fid,'W',h5.type_id,h5.space_id,h5.dcpl);

H5S.close(h5.space_id);
H5T.close(h5.type_id);
H5D.close(h5.dset_id);
H5F.close(h5.fid);

h5disp(hdf5File);