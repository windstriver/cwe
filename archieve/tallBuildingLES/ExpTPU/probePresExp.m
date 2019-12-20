% Convert the location of measured points from TPU experiment to
% the pressure probes used in OpenFOAM
% The building bottom center is located at the origin of the
% coordinate system of the computational domain

load('T215_4_000_1.mat');
% building height
H = Building_height;
% building breadth
B = Building_breadth;
% building depth
D = Building_depth;

% location of measured points
expPtLoc = Location_of_measured_points';
% number of measured points
nPt = size(expPtLoc, 1);

% Pressure probes coordinates
probeLoc = zeros(nPt, 5);
for i = 1:nPt
    x = expPtLoc(i,1);
    y = expPtLoc(i,2);
    id = expPtLoc(i,3);
    idSurface = expPtLoc(i,4);
    if (idSurface==1)    % surface 1: windward
        probeLoc(i,1) = -D/2;
        probeLoc(i,2) = B/2 - x;
        probeLoc(i,3) = y;
    end
    if (idSurface==2)    % surface 2: R-sideward
        probeLoc(i,1) = x - B - D/2;
        probeLoc(i,2) = -B/2;
        probeLoc(i,3) = y;
    end
    if (idSurface==3)    % surface 3: leeward
        probeLoc(i,1) = D/2;
        probeLoc(i,2) = x - (B+D) - B/2;
        probeLoc(i,3) = y;
    end
    if (idSurface==4)    % surface 4: L-sideward
        probeLoc(i,1) = D/2 - (x - (2*B+D));
        probeLoc(i,2) = B/2;
        probeLoc(i,3) = y;
    end
    probeLoc(i,4) = id;
    probeLoc(i,5) = idSurface;
end

% Write to CSV files
csvwrite('Location_of_CFD_pressure_probes.csv', probeLoc);
csvwrite('Location_of_measured_points.csv', expPtLoc);

% Location of pressure probes used in OpenFOAM
fileID = fopen('probePresBuildingExp', 'w');
for i = 1:nPt
    fprintf(fileID, '(% 6.4f  % 6.4f  % 6.4f)\n', ...
        probeLoc(i,1), probeLoc(i,2), probeLoc(i,3));
end
fclose(fileID);
