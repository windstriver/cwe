% Generate probe coordinates to sample approaching flow
par.B = 0.628;
par.H = 0.060;

fileID = fopen('probes', 'w');

for x = (-8:1:10)*par.B
    for y = (10:14)*par.H
        fprintf(fileID, '(% 8.4f % 8.4f % 8.4f)\n', x, y, 0);
    end
end

fclose(fileID);
    
