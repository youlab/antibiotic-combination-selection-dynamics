%% mantis_generate_2122020.m

% designs conditions and scrambles them for plate layout. writes protocol
% to protocolfile (specified each time).

prompt = "This script will overwrite your protocol file. Do you want to proceed? Y/N: ";
proceed = input(prompt, 's');
if ~strcmp(proceed, 'Y')
    return
end

% scramble rng
rng('shuffle');

%% Plate Information

% plate dimensions
height = 16;
width = 24;
nblanks = 96;
ninternalblanks = 20;
nwells = height * width - nblanks;
nreplicates = 4; % technical replicates
nconditions = nwells / nreplicates;
nisolates = 8;
nbioreplicates = 3;
% conditions = 72 = 8 isolates * 3 biological replicates * 3 conditions

%% create matrices

% read in condition details
protocolfile = '2-14-2020 3194 3398 Protocol.xlsx';
conditionmap = readmatrix(protocolfile, 'Sheet', 2, 'Range', 'A:E');

% condition % isolate %br %antibiotic concentration %inhibitor concentration
% conditionmap = [(1:nconditions).', zeros(nconditions, 4)];
% conditionmap(:, 2) = repmat(1:nisolates, 1, nconditions / nisolates);

% make and shuffle condition matrix
condition_replicates = repmat(1:nconditions, 1, nreplicates);
% add blanks and controls
conditions_withblanks = cat(2, condition_replicates, zeros(1, ninternalblanks));
conditions_shuffled = conditions_withblanks(randperm(size(conditions_withblanks, 2)));
conditions = reshape(conditions_shuffled, height-2, width-2);
fullconditions = zeros(height, width);
fullconditions(2:height-1, 2:width-1) = conditions;

writematrix(fullconditions, protocolfile, 'Sheet', 2, 'Range', 'G1:AD16');

%% make maps for mantis

% will end up with 24 isolate maps, a medium map, an antibiotic map, and an
% inhibitor map

isolates = zeros(nisolates, nbioreplicates, height, width);
antibiotic = zeros(height, width);
inhibitor = zeros(height, width);
medium = 100 * ones(height, width);

for i = 1:height
    for j = 1:width
        if fullconditions(i, j) > 0 % everything but blanks
        mapping = conditionmap(fullconditions(i, j), :);
        isolates(mapping(2), mapping(3), i, j) = 10; % cells
        antibiotic(i, j) = mapping(4);
        inhibitor(i, j) = mapping(5);
        end
    end
end

cells = squeeze(sum(isolates, 1));
cells = squeeze(sum(cells, 1));

medium = zeros(height, width);
medium = 100 - (cells + antibiotic + inhibitor);

writematrix(medium, protocolfile, 'Sheet', 3);
writematrix(medium, protocolfile, 'Sheet', 1, 'Range', 'BC8:BZ23');
writematrix(antibiotic, protocolfile, 'Sheet', 4);
writematrix(antibiotic, protocolfile, 'Sheet', 1, 'Range', 'A8:X23');
writematrix(inhibitor, protocolfile, 'Sheet', 5);
writematrix(inhibitor, protocolfile, 'Sheet', 1, 'Range', 'A26:X41');
%% write isolates out
counter = 6;
for i = 1:nisolates
    for j = 1:nbioreplicates 
        rangestart = 18*((nbioreplicates*(i-1)+j)-1)+44;
        rangeend = rangestart + 15;
        
        range = strcat('A', num2str(rangestart), ':X', num2str(rangeend));

        writematrix(squeeze(isolates(i, j, :, :)), protocolfile, 'Sheet', counter);
        counter = counter + 1;
        
        writematrix(squeeze(isolates(i, j, :, :)), protocolfile, 'Sheet', 1, 'Range', range);
    end
end

%%

% total volumes

% volume_antibiotic_low = sum(sum(antibiotic_low));
% volume_antibiotic_high = sum(sum(antibiotic_high));
% volume_inhibitor_low = sum(sum(inhibitor_low));
% volume_inhibitor_high = sum(sum(inhibitor_high));
% volume_medium_low = sum(sum(medium_low));
% volume_medium_high = sum(sum(medium_high));
% volume_culture = sum(sum(culture));
