%% mantis_recover.m: the save if you accidentally somehow overwrite your protocolfile but can recover the protocol via the Mantis dispense list

protocolfile = '1-28-2020 2165 2235 Protocol.xlsx';
conditionmap = readmatrix(protocolfile, 'Sheet', 2, 'Range', 'A:E');

height = 16;
width = 24;
nblanks = 96;
ninternalblanks = 20;
nwells = height * width - nblanks;
nreplicates = 4; % technical replicates
nconditions = nwells / nreplicates;
nisolates = 8;
nbioreplicates = 3;

isolates = zeros(nisolates, nbioreplicates, height, width);

% read everything back in
medium = readmatrix(protocolfile, 'Sheet', 1, 'Range', 'BC8:BZ23');
antibiotic = readmatrix(protocolfile, 'Sheet', 1, 'Range', 'A8:X23');
inhibitor = readmatrix(protocolfile, 'Sheet', 1, 'Range', 'A26:X41');

for i = 1:nisolates
    for j = 1:nbioreplicates 
        rangestart = 18*((nbioreplicates*(i-1)+j)-1)+44;
        rangeend = rangestart + 15;
        
        range = strcat('A', num2str(rangestart), ':X', num2str(rangeend));

        isolates(i, j, :, :) = readmatrix(protocolfile, 'Sheet', 1, 'Range', range);
    end
end

fullconditions = zeros(height, width);

for i = 1:height
    for j = 1:width
        key = zeros(1, 4); % mapping has 4 components: isolate, biorep, antibiotic, inhibitor
        key(4) = inhibitor(i, j);
        key(3) = antibiotic(i, j);
        
        for m = 1:nisolates
            for n = 1:nbioreplicates
                if isolates(m, n, i, j) == 10
                    key(1) = m;
                    key(2) = n;
                end
            end
        end
        
        % decode key
        for z = 1:nconditions
            if key == conditionmap(z, 2:5)
                fullconditions(i, j) = conditionmap(z, 1);
            end
        end
    end
end