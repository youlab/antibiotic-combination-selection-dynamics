%% mantis_generate_1312021.m

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
nwells = 384 - (2*height + 2*width - 4);
wellvol = 100;
cellvol = 20;
avol = (wellvol-cellvol)/2;
ivol = avol;

% experimental conditions
antibiotics = [0 0.5 1 2 4 8 16 32 64 128];
inhibitors = antibiotics;

nconditions = length(antibiotics)*length(inhibitors);

% calculate replicates and blanks
nreplicates = floor(nwells / nconditions);
ntotalblanks = 384 - (nconditions*nreplicates);
ninternalblanks = ntotalblanks - (2*height + 2*width - 4);


%% create matrices

% read in condition details
protocolfile = '6-23-2021 Protocol.xlsx';
conditionmap = readmatrix(protocolfile, 'Sheet', 2, 'Range', 'A2:C101'); % changed because first row in conditionmap reading NaN (row titles)
% condition % strain % plasmid % medium

% make and shuffle condition matrix
condition_replicates = repmat(1:nconditions, 1, nreplicates);
% add blanks and controls
conditions_withblanks = cat(2, condition_replicates, zeros(1, ninternalblanks));
conditions_shuffled = conditions_withblanks(randperm(size(conditions_withblanks, 2)));
conditions = reshape(conditions_shuffled, height-2, width-2);
fullconditions = zeros(height, width);
fullconditions(2:height-1, 2:width-1) = conditions;

writematrix(fullconditions, protocolfile, 'Sheet', 2, 'Range', 'E1:AB16');

%% make maps for mantis

amaps = zeros(length(antibiotics), height, width);
imaps = amaps;
cells = zeros(height, width);
blanks = zeros(height, width);

for i = 1:height
    for j = 1:width
        if fullconditions(i, j) > 0 % everything but blanks
            mapping = conditionmap(fullconditions(i, j), :);
            amaps(mapping(2), i, j) = avol;
            imaps(mapping(3), i, j) = ivol;
            cells(i, j) = cellvol;
        elseif fullconditions(i, j) == 0 % add blanks
            blanks(i, j) = wellvol;
        end
    end
end

writematrix(blanks, protocolfile, 'Sheet', 3);
writematrix(blanks, protocolfile, 'Sheet', 1, 'Range', 'A9:X24');
writematrix(cells, protocolfile, 'Sheet', 4);
writematrix(cells, protocolfile, 'Sheet', 1, 'Range', 'A408:X423');

%% write out maps

counter = 5;

for i = 1:length(antibiotics)
    rangestart = 19*i + 9;
    rangeend = rangestart + 15;
    
    range = strcat('A', num2str(rangestart), ':X', num2str(rangeend));

    writematrix(squeeze(amaps(i, :, :)), protocolfile, 'Sheet', counter);
    counter = counter + 1;
    
    writematrix(squeeze(amaps(i, :, :)), protocolfile, 'Sheet', 1, 'Range', range);

end

for i = 1:length(inhibitors)
    rangestart = 19*i+199;
    rangeend = rangestart + 15;

    range = strcat('A', num2str(rangestart), ':X', num2str(rangeend));

    writematrix(squeeze(imaps(i, :, :)), protocolfile, 'Sheet', counter);
    counter = counter + 1;

    writematrix(squeeze(imaps(i, :, :)), protocolfile, 'Sheet', 1, 'Range', range);

end