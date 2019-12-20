% Set the pressure probes before the building to monitor the propagation
% of the inflow turbulence

load('T215_4_000_1.mat');
% building height
H = Building_height;
% building breadth
B = Building_breadth;
% building depth
D = Building_depth;

xGrid = -5:0.2:0;
yGrid = -2.2:0.2:2.2;
zGrid = 0:0.2:1;

noGrid = length(xGrid) * length(yGrid) * length(zGrid);

probeLoc = zeros(noGrid, 7);
ip = 1;
for ix = 1:length(xGrid)
    for iy = 1:length(yGrid)
        for iz = 1:length(zGrid)
            probeLoc(ip,1) = xGrid(ix) - D/2;
            probeLoc(ip,2) = yGrid(iy);
            probeLoc(ip,3) = zGrid(iz);
            probeLoc(ip,4) = ip;
            probeLoc(ip,5) = ix;
            probeLoc(ip,6) = iy;
            probeLoc(ip,7) = iz;
            
            ip = ip + 1;
        end
    end
end

csvwrite('Location_of_CFD_velocity_probes.csv', probeLoc);

% Location of pressure probes used in OpenFOAM
fileID = fopen('probeUUpstream', 'w');
for i = 1:noGrid
    fprintf(fileID, '(% 6.4f  % 6.4f  % 6.4f)\n', ...
        probeLoc(i,1), probeLoc(i,2), probeLoc(i,3));
end
fclose(fileID);
