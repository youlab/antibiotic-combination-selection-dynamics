clear
close all

%% this file only supports a single endpoint and a single replicate

datafile = '3-19-2022 SUL 3 HM areas.xlsx';

data = readmatrix(datafile, 'Sheet', 1, 'Range', 'A:C');
info = readmatrix(datafile, 'Sheet', 1, 'Range', 'D:H');

%%

initialdata = data(info(:, 1) == 0, :);
data = data(~(info(:, 1) == 0), :);
info = info(~(info(:, 1) == 0), :);
%% split data out by condition
abxconcentrations = unique(info(~isnan(info(:, 2)), 2));
inhconcentrations = unique(info(~isnan(info(:, 3)), 3));
nimages = 3;
timepoints = unique(info(:, 1));

%timepoint [A] [I] image sensitive/resistant/fraction
datasplit = nan(length(abxconcentrations), length(inhconcentrations), nimages, size(data, 2));

for i = 1:size(data, 1)
    if ~isnan(info(i, 2))
        datasplit(abxconcentrations==info(i, 2), inhconcentrations==info(i, 3), info(i, 5), :) = data(i, :);
    end
end

imagessum = squeeze(sum(datasplit(:, :, :, 1:2), 3, 'omitnan'));
fractions = imagessum(:, :, 2) ./ imagessum(:, :, 1);

initialsum = squeeze(sum(initialdata(:, 1:2), 1, 'omitnan'));
initialfraction = initialsum(2) ./ initialsum(1);

alltimepoints = cat(1, [0], timepoints);


%% average fractions long

fractionslong = zeros(length(timepoints)+1*length(abxconcentrations)*length(inhconcentrations), 4);

for i = 1:length(alltimepoints)
    for j = 1:length(abxconcentrations)
        for k = 1:length(inhconcentrations)
            index = k + length(inhconcentrations)*(j-1) + (length(inhconcentrations) * length(abxconcentrations))*(i-1);
            if i == 1
                thisfraction = initialfraction;
            else
                thisfraction = fractions(j, k);
            end
            fractionslong(index, :) = [alltimepoints(i) abxconcentrations(j) inhconcentrations(k) thisfraction];
        end
    end
end

%save('SUL 3 HM microscopy.mat', 'fractionslong')
